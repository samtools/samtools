/*  test/fuzz/fuzz_common.h -- Shared code for samtools subcommand fuzzers

    Copyright (C) 2026 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*
Shared harness for fuzzing samtools subcommands through their main
functions.  Include this once, from the fuzzer's single source file.

A fuzzer supplies an argument builder and calls fuzz_one_input():

    #include "fuzz_common.h"
    int bam_flagstat(int argc, char *argv[]);

    static void args(fuzz_args *a, const uint8_t *c, char *in) {
        ARG("flagstat");
        CHOOSE("-O", c[0] & 3, NULL, "default", "json", "tsv");
        ENDOPTS();
        ARG(in);
    }

    int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
        return fuzz_one_input(data, size, 1, FUZZ_PARSE, args, bam_flagstat);
    }

Inputs are laid out as [ctrl_bytes control bytes][alignment data].  The
data is passed to samtools as a "data:" URL, so no input file is written.

Commands needing an index use fuzz_one_file_input() instead, which writes
the input to a file, indexes it, and generates regions; see the "File-backed
inputs" section below.

Build every fuzzer with -Dexit=fuzz_exit so exit() calls unwind back here.

Run as, e.g.:
ASAN_OPTIONS=allow_addr2line=1:abort_on_error=1:detect_leaks=0 ./fuzz_X -rss_limit_mb=2560 -timeout=10 -close_fd_mask=1

-close_fd_mask=1 discards stdout, where outputs go.  Set FUZZ_QUIET=1 in the
environment to stop each command line being echoed to stderr.

Reproducing a crash outside the harness:

    FUZZ_EXTRACT=repro ./fuzz_X crash-file
    cd repro && SAMTOOLS=/path/to/samtools ./repro.sh

Everything a run uses comes from the fuzz input itself, apart from a few
fixed runtime files (ref.fa etc.).  With FUZZ_EXTRACT=DIR the harness
writes every derived file into DIR under stable names ("fz.*") instead of
the temporary "_fz.<pid>.*" ones, keeps them, copies in any runtime files
the commands reference, and writes DIR/repro.sh containing the exact
command line(s), with paths relative to DIR.  The raw payloads (fz.input,
and fz.extra where used) are written first, so they exist even if the crash
is in the harness's own preparation.  The commands are then run as normal,
so the crash also reproduces in the harness; set FUZZ_EXTRACT_ONLY=1 to
just extract.  Intended for single-input runs.

Crashes that don't reproduce from a single input:

The harness calls each subcommand's main() repeatedly in one process, so
a crash can depend on state left behind by earlier inputs (static
variables, abandoned allocations after an exit() longjmp, ...).  libFuzzer
does not record the sequence, so set FUZZ_HISTORY=DIR while fuzzing: each
input is saved, before it runs, as DIR/hist-<pid>-<seq>, keeping the most
recent FUZZ_HISTORY_LEN (default 100) per process.  After a crash, replay
that process's history in order, in one process:

    ./fuzz_X -detect_leaks=0 $(ls HISTDIR/hist-<pid>-* | sort)

then bisect by dropping files from the front of the list.  A crash that
needs a particular earlier input usually points at state that a real
single-shot samtools never sees, but which may still be a genuine bug if
it lives in htslib, as library users call it repeatedly.

Also pass -detect_leaks=0 when reproducing: otherwise libFuzzer reruns an
input that doesn't balance its mallocs and frees, to look for a leak, so
everything appears twice.

Argument rules, so that fuzzing can never cause damage:

  1. Every argument is a string literal, a fixed name built from
     fuzz_prefix, or explicitly marked as data-derived: the input URL, or
     strings added with ARG_DATA() / ARG_OPTDATA() (e.g. region strings
     built from the fuzzed header).
  2. Options end with ENDOPTS() ("--") before any positional argument.
     A data-derived string may only appear after "--" (ARG_DATA), or as
     the value of an option that takes a required argument (ARG_OPTDATA,
     e.g. -r REGION), where getopt consumes it as the value regardless of
     its content.  The driver checks this on every run and aborts on any
     violation.
  3. No output path is fuzzable.  Outputs go to stdout, "-", /dev/null, or
     a fixed name under fuzz_prefix.  Options that build filenames from
     the data (e.g. stats --split, which uses tag values) are not used.
     Data-derived *content* may be written to fixed names (e.g. a BED
     file of generated regions).

Files: anything a fuzzer has to write (sort/markdup temporaries, stats
--split files) uses the prefix fuzz_prefix ("_fz.<pid>") in the current
directory, and everything matching that prefix is removed before and after
each run.  Run from a scratch directory.  Output that would otherwise go to
a named file is sent to /dev/null.

Runtime input files (ref.fa, regions.bed, ...) are fixed names in the
current directory; see make_fuzz_data.py.
*/

#ifndef FUZZ_COMMON_H
#define FUZZ_COMMON_H

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>
#include <getopt.h>
#include <glob.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <htslib/sam.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <inttypes.h>

// Taken from bamtk.c which we exclude from our build line
const char *samtools_version(void) {
    return "fuzzer";
}

/* ------------------------------------------------------------------------
 * Input encoding and logging
 */

// Returns "data:;base64,<base64-data>" as a malloced buffer
static char *fuzz_base64(const uint8_t *data, size_t len) {
    static const char table[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    char *output = malloc(4 * ((len + 2) / 3) + 20);
    if (!output)
        return NULL;

    strcpy(output, "data:;base64,");
    size_t i = 0, j = strlen(output);

    // Padding is decided by the bytes remaining *before* each group.
    while (i < len) {
        size_t rem = len - i;
        unsigned int a =             data[i++];
        unsigned int b = (rem > 1) ? data[i++] : 0;
        unsigned int c = (rem > 2) ? data[i++] : 0;
        unsigned int triple = (a << 16) | (b << 8) | c;

        output[j++] = table[(triple >> 18) & 0x3F];
        output[j++] = table[(triple >> 12) & 0x3F];
        output[j++] = (rem > 1) ? table[(triple >> 6) & 0x3F] : '=';
        output[j++] = (rem > 2) ? table[triple & 0x3F]        : '=';
    }

    output[j] = '\0';
    return output;
}

// Writes one shell word, quoted if necessary
static void fuzz_quote(FILE *fp, const char *a) {
    static const char safe[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz0123456789_-.,:/=+@%";
    if (*a && strspn(a, safe) == strlen(a)) {
        fputs(a, fp);
        return;
    }
    fputc('\'', fp);
    for (; *a; a++)
        if (*a == '\'') fputs("'\\''", fp);
        else             fputc(*a, fp);
    fputc('\'', fp);
}

// Prints the command with shell quoting, so it can be pasted to reproduce.
static void fuzz_dump(int argc, char **argv) {
    fprintf(stderr, "Calling:");
    for (int i = 0; i < argc; i++) {
        fputc(' ', stderr);
        fuzz_quote(stderr, argv[i]);
    }
    fputc('\n', stderr);
}

/* ------------------------------------------------------------------------
 * exit() trapping and initialisation
 */

static jmp_buf fuzz_jmp;
void fuzz_exit(int status) {
    fprintf(stderr, "Caught exit of status %d\n", status);
    longjmp(fuzz_jmp, status + 1);
}

// Prefix for any files a fuzzer has to write, removed after each run
// (unless extracting)
static char fuzz_prefix[300];
static int fuzz_quiet;

// FUZZ_EXTRACT directory, or "" when not extracting.  When extracting, the
// harness chdir()s there on the first run (libFuzzer reads its input files
// after initialisation, so not before), so every path is the same for the
// harness and for repro.sh.  Runtime files are copied in from fuzz_orig_cwd.
static char fuzz_extract_dir[256];
static char fuzz_orig_cwd[4096];
static int fuzz_extract_only, fuzz_extract_entered;

// FUZZ_HISTORY: an absolute directory, or "" when not recording
static char fuzz_hist_dir[4096];
static unsigned long long fuzz_hist_seq;
static unsigned long long fuzz_hist_len = 100;

int LLVMFuzzerInitialize(int *argc, char ***argv) {
    (void) argc; (void) argv;
    // CRAM input would otherwise fetch references from the EBI server
    // and cache them under ~/.cache/hts-ref.  Make lookups local-only
    // and doomed to fail, so we never touch the network or write files.
    setenv("REF_PATH", "_fuzz_no_ref_path_", 1);
    unsetenv("REF_CACHE");
    snprintf(fuzz_prefix, sizeof(fuzz_prefix), "_fz.%d", (int) getpid());
    fuzz_quiet = getenv("FUZZ_QUIET") != NULL;

    const char *ex = getenv("FUZZ_EXTRACT");
    if (ex && *ex) {
        if (strlen(ex) >= sizeof(fuzz_extract_dir) - 8) {
            fprintf(stderr, "FUZZ_EXTRACT path too long\n");
            abort();
        }
        if (mkdir(ex, 0777) < 0 && errno != EEXIST) {
            perror(ex);
            abort();
        }
        snprintf(fuzz_extract_dir, sizeof(fuzz_extract_dir), "%s", ex);
        snprintf(fuzz_prefix, sizeof(fuzz_prefix), "fz");
        if (!getcwd(fuzz_orig_cwd, sizeof(fuzz_orig_cwd))) {
            perror("getcwd");
            abort();
        }
        fuzz_extract_only = getenv("FUZZ_EXTRACT_ONLY") != NULL;
    }

    const char *h = getenv("FUZZ_HISTORY");
    if (h && *h) {
        if (mkdir(h, 0777) < 0 && errno != EEXIST) {
            perror(h);
            abort();
        }
        char cwd[2048];
        if (h[0] != '/' && !getcwd(cwd, sizeof(cwd))) {
            perror("getcwd");
            abort();
        }
        if (h[0] == '/')
            snprintf(fuzz_hist_dir, sizeof(fuzz_hist_dir), "%s", h);
        else
            snprintf(fuzz_hist_dir, sizeof(fuzz_hist_dir), "%s/%s", cwd, h);
        const char *n = getenv("FUZZ_HISTORY_LEN");
        if (n && atoll(n) > 0)
            fuzz_hist_len = atoll(n);
    }
    return 0;
}

static void fuzz_cleanup(void) {
    char pattern[320];
    glob_t g;
    snprintf(pattern, sizeof(pattern), "%s*", fuzz_prefix);
    if (glob(pattern, 0, NULL, &g) == 0) {
        for (size_t i = 0; i < g.gl_pathc; i++)
            unlink(g.gl_pathv[i]);
        globfree(&g);
    }
}

/* ------------------------------------------------------------------------
 * Input pre-checks
 */

#define FUZZ_PARSE        0  // reject inputs that fail to parse
#define FUZZ_SORTED       1  // ... or are not coordinate sorted
#define FUZZ_SHORT_REFS   2  // ... or have @SQ lengths over 1000 (for -a)
#define FUZZ_ANYFMT       4  // accept anything: inputs that don't parse as
                             // alignments (FASTA, FASTQ, damaged files) are
                             // passed through as they are

// Reject inputs that don't parse, as we don't need to fuzz that (htslib
// does already), and our tests are slow so reserve them for correct data.
// Returns 0 if usable, -1 if not.
static int fuzz_check_input(const uint8_t *data, size_t size, int flags) {
    uint8_t *copy = malloc(size);
    if (!copy)
        abort();
    memcpy(copy, data, size);

    hFILE *memfile = hopen("mem:", "rb:", copy, size);
    if (!memfile) {
        free(copy);
        return -1;
    }
    htsFile *in = hts_hopen(memfile, "data", "rb");
    if (!in) {
        if (hclose(memfile) != 0)
            abort();
        return -1;
    }

    sam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
        hts_close(in);
        return -1;
    }

    int ret = 0;
    if (flags & FUZZ_SHORT_REFS) {
        for (int i = 0; i < sam_hdr_nref(hdr); i++) {
            if (sam_hdr_tid2len(hdr, i) > 1000) {
                ret = -1;
                break;
            }
        }
    }

    bam1_t *b = bam_init1();
    if (!b)
        abort();
    int r = 0;
    int64_t last_tid = -1;
    hts_pos_t last_pos = -1;
    while (ret == 0 && (r = sam_read1(in, hdr, b)) >= 0) {
        if (flags & FUZZ_SORTED) {
            // Unplaced reads (tid -1) sort last
            int64_t tid = b->core.tid < 0 ? INT64_MAX : b->core.tid;
            if (tid < last_tid || (tid == last_tid && b->core.pos < last_pos))
                ret = -1;
            last_tid = tid;
            last_pos = b->core.pos;
        }
    }
    if (r < -1)
        ret = -1; // parse error

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    hts_close(in);
    return ret;
}

/* ------------------------------------------------------------------------
 * Argument building
 */

#define FUZZ_MAX_ARGS 200

// How an argument may relate to the fuzz data
#define FUZZ_FIXED    0  // literal or fixed name
#define FUZZ_DATA     1  // data-derived; must follow "--"
#define FUZZ_OPTDATA  2  // data-derived value of the preceding option

typedef struct {
    char *argv[FUZZ_MAX_ARGS];
    unsigned char kind[FUZZ_MAX_ARGS];
    int argc;
} fuzz_args;

static void fuzz_arg_kind(fuzz_args *a, const char *s, int kind) {
    if (a->argc >= FUZZ_MAX_ARGS - 1)
        abort(); // harness bug: raise FUZZ_MAX_ARGS
    a->kind[a->argc] = kind;
    a->argv[a->argc++] = (char *) s;
}

static void fuzz_arg(fuzz_args *a, const char *s) {
    fuzz_arg_kind(a, s, FUZZ_FIXED);
}

// Enforces rule 2.  "input" (may be NULL) is also treated as data-derived.
static void fuzz_check_args(const fuzz_args *a, const char *input) {
    int i, endopts = -1;
    for (i = 1; i < a->argc; i++) {
        int kind = a->kind[i];
        if (input && a->argv[i] == input && kind == FUZZ_FIXED)
            kind = FUZZ_DATA;
        if (endopts < 0 && kind == FUZZ_FIXED && strcmp(a->argv[i], "--") == 0)
            endopts = i;
        if (kind == FUZZ_DATA && endopts < 0) {
            fprintf(stderr, "Harness bug: data-derived argument %d "
                    "before \"--\"\n", i);
            abort();
        }
        if (kind == FUZZ_OPTDATA &&
            (endopts >= 0 || a->kind[i-1] != FUZZ_FIXED
             || a->argv[i-1][0] != '-')) {
            fprintf(stderr, "Harness bug: data-derived argument %d is not "
                    "an option value\n", i);
            abort();
        }
    }
    if (endopts < 0) {
        fprintf(stderr, "Harness bug: no \"--\"\n");
        abort();
    }
}

// Picks vals[sel % n]; a NULL value omits the option entirely.  With a
// NULL opt, the value alone is added (for a choice of flags).
static void fuzz_choose(fuzz_args *a, const char *opt, unsigned sel,
                        const char *const *vals, size_t n) {
    const char *v = vals[sel % n];
    if (!v)
        return;
    if (opt)
        fuzz_arg(a, opt);
    fuzz_arg(a, v);
}

// These assume the builder's fuzz_args pointer is called "a".
#define ARG(s)       fuzz_arg(a, (s))
#define ENDOPTS()    fuzz_arg(a, "--")
#define ARG_DATA(s)  fuzz_arg_kind(a, (s), FUZZ_DATA)
#define ARG_OPTDATA(o, v) \
    do { fuzz_arg(a, (o)); fuzz_arg_kind(a, (v), FUZZ_OPTDATA); } while (0)
#define ARG2(o, v)   do { fuzz_arg(a, (o)); fuzz_arg(a, (v)); } while (0)
#define CHOOSE(opt, sel, ...)                                              \
    fuzz_choose(a, (opt), (sel), (const char *const []){__VA_ARGS__},       \
                sizeof((const char *const []){__VA_ARGS__}) / sizeof(char *))

/* ------------------------------------------------------------------------
 * Extraction (FUZZ_EXTRACT)
 */

static void fuzz_write_fixed(const char *fn, const void *data, size_t len) {
    FILE *fp = fopen(fn, "wb");
    if (!fp || fwrite(data, 1, len, fp) != len || fclose(fp) != 0) {
        perror(fn); // infrastructure failure, not a samtools bug
        abort();
    }
}

// Moves into the extraction directory on the first run
static void fuzz_extract_enter(void) {
    if (!*fuzz_extract_dir || fuzz_extract_entered)
        return;
    if (chdir(fuzz_extract_dir) < 0) {
        perror(fuzz_extract_dir);
        abort();
    }
    fuzz_extract_entered = 1;
}

// Saves a raw payload into the extraction directory, if extracting
static void fuzz_extract_payload(const char *name, const void *data,
                                 size_t len) {
    if (*fuzz_extract_dir)
        fuzz_write_fixed(name, data, len);
}

// Copies a runtime file (a bare name in the original directory) into the
// extraction directory, along with its .fai/.gzi companions
static void fuzz_extract_runtime(const char *name) {
    static const char *sfx[3] = {"", ".fai", ".gzi"};
    struct stat st;
    if (!*name || strchr(name, '/') || strcmp(name, "-") == 0
        || strncmp(name, fuzz_prefix, strlen(fuzz_prefix)) == 0)
        return;
    for (int i = 0; i < 3; i++) {
        char src[4400], dst[300];
        snprintf(src, sizeof(src), "%s/%s%s", fuzz_orig_cwd, name, sfx[i]);
        snprintf(dst, sizeof(dst), "%s%s", name, sfx[i]);
        if (stat(src, &st) < 0 || !S_ISREG(st.st_mode))
            continue;
        FILE *in = fopen(src, "rb"), *out = fopen(dst, "wb");
        if (!in || !out) {
            perror(!in ? src : dst);
            abort();
        }
        char buf[8192];
        size_t n;
        while ((n = fread(buf, 1, sizeof(buf), in)) > 0)
            if (fwrite(buf, 1, n, out) != n)
                abort();
        fclose(in);
        if (fclose(out) != 0)
            abort();
    }
}

// Writes an argument to the repro script.  "input" (if non-NULL) is the
// data: URL of fz.input, which the script regenerates from that file.
static void fuzz_repro_arg(FILE *fp, const char *arg, const char *input) {
    if (input && arg == input) {
        fputs("\"data:;base64,$(base64 < fz.input | tr -d '\\n')\"", fp);
        return;
    }
    fuzz_quote(fp, arg);
    fuzz_extract_runtime(arg);
}

// Writes DIR/repro.sh for these command lines (replacing any earlier one)
static void fuzz_write_repro(const fuzz_args *cmds, int ncmds,
                             const char *input) {
    const char *fn = "repro.sh";
    FILE *fp = fopen(fn, "w");
    if (!fp) {
        perror(fn);
        abort();
    }
    fprintf(fp, "#!/bin/sh\n"
            "# Generated by the samtools fuzz harness from a fuzz input.\n"
            "# Set SAMTOOLS to choose the binary, e.g. an ASAN build.\n"
            "cd \"$(dirname \"$0\")\" || exit 1\n"
            "SAMTOOLS=${SAMTOOLS:-samtools}\n");
    for (int c = 0; c < ncmds; c++) {
        // Commands other than samtools (e.g. an "rm" between steps) are
        // written as they are
        if (strcmp(cmds[c].argv[0], "samtools") == 0)
            fputs("\"$SAMTOOLS\"", fp);
        else
            fuzz_quote(fp, cmds[c].argv[0]);
        for (int i = 1; i < cmds[c].argc; i++) {
            fputc(' ', fp);
            fuzz_repro_arg(fp, cmds[c].argv[i], input);
        }
        fputc('\n', fp);
    }
    if (fclose(fp) != 0)
        abort();
    chmod(fn, 0755);
    fprintf(stderr, "Extracted to %s/; run %s/repro.sh\n",
            fuzz_extract_dir, fuzz_extract_dir);
}

/* ------------------------------------------------------------------------
 * Run hygiene: input history, and the standard file descriptors
 */

// Saves the whole input (control bytes included) before it runs, keeping
// the most recent fuzz_hist_len per process
static void fuzz_record_history(const uint8_t *data, size_t size) {
    if (!*fuzz_hist_dir)
        return;
    char fn[4200];
    int pid = (int) getpid();
    snprintf(fn, sizeof(fn), "%s/hist-%d-%012llu", fuzz_hist_dir, pid,
             fuzz_hist_seq);
    fuzz_write_fixed(fn, data, size);
    if (fuzz_hist_seq >= fuzz_hist_len) {
        snprintf(fn, sizeof(fn), "%s/hist-%d-%012llu", fuzz_hist_dir, pid,
                 fuzz_hist_seq - fuzz_hist_len);
        unlink(fn);
    }
    fuzz_hist_seq++;
}

// A subcommand that closes stdin/stdout/stderr (e.g. by closing a "-"
// output) would otherwise leave the next run writing to whatever file
// later reuses that descriptor.  Snapshot 0-2 on the first run (after
// libFuzzer's -close_fd_mask redirection), and restore them before and
// after every command.
static int fuzz_std_fd[3] = {-1, -1, -1};
static struct stat fuzz_std_st[3];

static void fuzz_std_fds(void) {
    fflush(stdout);
    fflush(stderr);
    for (int i = 0; i < 3; i++) {
        struct stat st;
        if (fuzz_std_fd[i] < 0) {
            // First call: take the snapshot, well clear of low fds
            if ((fuzz_std_fd[i] = fcntl(i, F_DUPFD_CLOEXEC, 100)) >= 0)
                fstat(fuzz_std_fd[i], &fuzz_std_st[i]);
            continue;
        }
        if (fstat(i, &st) < 0 || st.st_dev != fuzz_std_st[i].st_dev
            || st.st_ino != fuzz_std_st[i].st_ino)
            dup2(fuzz_std_fd[i], i);
    }
    clearerr(stdout);
    clearerr(stderr);
}

static void fuzz_begin_run(const uint8_t *data, size_t size) {
    fuzz_record_history(data, size);
    fuzz_std_fds();
}

/* ------------------------------------------------------------------------
 * Writable argv
 *
 * A real main() gets writable argv strings, and some subcommands modify
 * their arguments in place (e.g. mpileup --output-extra via strtok_r).
 * String literals live in read-only memory, so every argument is copied
 * before the subcommand is called.  The copies are freed at the start of
 * the next call, so they're also reclaimed after an exit() longjmp.
 */

static char *fuzz_argv_copy[FUZZ_MAX_ARGS];

static void fuzz_free_argv(void) {
    for (int i = 0; i < FUZZ_MAX_ARGS && fuzz_argv_copy[i]; i++) {
        free(fuzz_argv_copy[i]);
        fuzz_argv_copy[i] = NULL;
    }
}

static char **fuzz_writable_argv(const fuzz_args *a) {
    fuzz_free_argv();
    for (int i = 0; i < a->argc; i++)
        if (!(fuzz_argv_copy[i] = strdup(a->argv[i])))
            abort();
    fuzz_argv_copy[a->argc] = NULL;
    return fuzz_argv_copy;
}

/* ------------------------------------------------------------------------
 * The driver
 */

typedef void (*fuzz_build_fn)(fuzz_args *a, const uint8_t *ctrl, char *input);
typedef int  (*fuzz_main_fn)(int argc, char *argv[]);

// Global so it can be freed on the next run if we longjmp'd out of this one
static char *fuzz_b64 = NULL;

static int fuzz_one_input(const uint8_t *data, size_t size, size_t ctrl_bytes,
                          int check_flags, fuzz_build_fn build,
                          fuzz_main_fn main_fn) {
    fuzz_begin_run(data, size);
    free(fuzz_b64);
    fuzz_b64 = NULL;
    fuzz_extract_enter();
    fuzz_cleanup();

    if (setjmp(fuzz_jmp) != 0) {
        // exit() was called; don't add it to the test corpus.
        // Files are tidied up at the start of the next run.
        return -1;
    }

    if (size < ctrl_bytes + 5)
        return 0;

    // Save the raw payload first, in case the crash is in the pre-check
    fuzz_extract_payload("fz.input", data + ctrl_bytes, size - ctrl_bytes);

    if (!(check_flags & FUZZ_ANYFMT) &&
        fuzz_check_input(data + ctrl_bytes, size - ctrl_bytes,
                         check_flags) < 0)
        return -1;

    fuzz_b64 = fuzz_base64(data + ctrl_bytes, size - ctrl_bytes);
    if (!fuzz_b64)
        abort();

    fuzz_args a = {{0}, {0}, 0};
    fuzz_arg(&a, "samtools");
    build(&a, data, fuzz_b64);

    // Enforce rule 2, and check the input was actually passed
    fuzz_check_args(&a, fuzz_b64);
    int i, n_input = 0;
    for (i = 1; i < a.argc; i++)
        n_input += a.argv[i] == fuzz_b64;
    if (n_input == 0) {
        fprintf(stderr, "Harness bug: input not passed\n");
        abort();
    }

    if (!fuzz_quiet)
        fuzz_dump(a.argc, a.argv);
    if (*fuzz_extract_dir) {
        fuzz_write_repro(&a, 1, fuzz_b64);
        if (fuzz_extract_only)
            return 0;
    }

#ifdef __GLIBC__
    optind = 0; // full reset, including glibc's internal short-option state
#else
    optind = 1;
#endif
    opterr = 0;

    char **argv = fuzz_writable_argv(&a);
    (void) main_fn(a.argc - 1, argv + 1);
    fuzz_free_argv();
    fuzz_std_fds();

    free(fuzz_b64);
    fuzz_b64 = NULL;
    if (!*fuzz_extract_dir)
        fuzz_cleanup();

    return 0;
}

/* ------------------------------------------------------------------------
 * File-backed inputs, for commands that need an index
 *
 * Inputs are laid out as:
 *   [1 file byte][tool_bytes tool bytes][FUZZ_MAX_REGIONS * FUZZ_REG_BYTES
 *   region bytes][alignment data]
 *
 * File byte:
 *   bits 0-1  how the data is written: BAM, CRAM (no_ref) or BGZF SAM, all
 *             re-encoded by htslib, or the raw input bytes (indexable only
 *             if the input itself is BAM/CRAM/BGZF SAM).  With FUZZ_ANYFMT,
 *             data that doesn't parse as alignments is always written raw,
 *             and gets regions with an absent reference name.
 *   bits 2-3  index: default name and type; CSI; alternate fixed name
 *             (f->idx_fn, for -X or "##idx##"); or no index at all
 *   bits 4-5  number of regions - 1
 *   bits 6-7  unused
 *
 * Region bytes (per region):
 *   r[0] bits 0-2  syntax: name, name:beg, name:beg-end, name:beg-,
 *                  {name}:beg-end, {name}, name:beg-<overflow>, special
 *   r[0] bits 3-7  reference selector (low bits)
 *   r[1], r[2]     start and end: near 0, near the end, past the end,
 *                  or proportional (top 2 bits pick which)
 *   r[3] bit  0    thousands separators
 *   r[3] bits 1-2  for "special": "*", ".", absent name, empty string
 *   r[3] bits 3-7  reference selector (high bits)
 *
 * With FUZZ_EXTRA in the flags, the region bytes are followed by a second
 * payload: [2-byte little-endian length L][L bytes], then the alignment
 * data.  The extra payload is written as given to f->extra_raw_fn (for use
 * as a header file, BED file, ...), and if it parses as alignment data it
 * is also written to f->extra_fn in the same form as f->fn, and made
 * available as a data: URL in f->extra_url, for a second input.
 *
 * The driver writes the data to f->fn, builds the index (with FUZZ_INDEX,
 * unless the fuzzer is testing index building itself), and then runs each
 * step in turn.  Builders pass f->in as the input, which is "fn##idx##idx_fn" when
 * the index has the alternate name, or f->fn and f->idx_path with -X.
 */

#define FUZZ_REG_BYTES   4
#define FUZZ_MAX_REGIONS 4

// fuzz_one_file_input() flags
#define FUZZ_INDEX  1   // build an index for the input(s)
#define FUZZ_EXTRA  2   // the input carries an extra payload

typedef struct {
    char fn[320];       // the alignment file
    char idx_fn[320];   // alternate index name
    char idx_path[340]; // the index actually built, for -X
    char idx_url[680];  // fn##idx##idx_fn
    char bed_fn[320];   // the generated regions as BED
    char list_fn[320];  // a file-of-filenames listing f->in
    char extra_raw_fn[320]; // FUZZ_EXTRA payload, as given
    char extra_fn[320]; // FUZZ_EXTRA payload re-encoded like fn
    char *extra_url;    // FUZZ_EXTRA payload as a data: URL (data-derived)
    int extra_aln;      // the extra payload parsed as alignment data
    const char *in;     // what to pass as the input: fn or idx_url
    int form;           // 0 BAM, 1 CRAM, 2 BGZF SAM, 3 raw
    int alt_idx;        // index uses idx_fn
    sam_hdr_t *hdr;
    kstring_t reg[FUZZ_MAX_REGIONS];
    int nreg;
} fuzz_file;

typedef void (*fuzz_file_build_fn)(fuzz_args *a, const uint8_t *ctrl,
                                   fuzz_file *f);
typedef struct {
    fuzz_file_build_fn build;
    fuzz_main_fn main_fn;   // NULL for a pseudo-command
    const char *pseudo;     // "rm" or "truncate", run by the harness
} fuzz_step;

// Runs a pseudo-command step, with the semantics of the shell command of
// the same name (as written to repro.sh):
//   rm -f -- FILE...
//   truncate -s [+-]N -- FILE...   (sign: relative to the current size)
static void fuzz_run_pseudo(const fuzz_args *a) {
    const char *cmd = a->argv[0], *size = NULL;
    int i;
    for (i = 1; i < a->argc && strcmp(a->argv[i], "--") != 0; i++)
        if (strcmp(a->argv[i], "-s") == 0 && i + 1 < a->argc)
            size = a->argv[++i];
    for (i++; i < a->argc; i++) {
        if (strcmp(cmd, "rm") == 0) {
            unlink(a->argv[i]);
        } else if (strcmp(cmd, "truncate") == 0 && size) {
            struct stat st;
            if (stat(a->argv[i], &st) < 0)
                continue;
            long long n = atoll(size);
            long long len = (*size == '-' || *size == '+') ? st.st_size + n : n;
            if (truncate(a->argv[i], len < 0 ? 0 : len) < 0)
                perror(a->argv[i]);
        } else {
            fprintf(stderr, "Harness bug: unknown pseudo-command %s\n", cmd);
            abort();
        }
    }
}

static fuzz_file fuzz_f;

// Parses alignment data and writes it to fn in the given form (0 BAM,
// 1 CRAM, 2 BGZF SAM, 3 raw bytes).  The header is returned in *hdr_out if
// non-NULL, otherwise freed.  Returns 0 on success.
static int fuzz_write_aln(const uint8_t *data, size_t size, int form,
                          const char *fn, sam_hdr_t **hdr_out) {
    uint8_t *copy = malloc(size);
    if (!copy)
        abort();
    memcpy(copy, data, size);
    hFILE *memfile = hopen("mem:", "rb:", copy, size);
    if (!memfile) {
        free(copy);
        return -1;
    }
    htsFile *in = hts_hopen(memfile, "data", "rb");
    if (!in) {
        hclose_abruptly(memfile);
        return -1;
    }
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
        hts_close(in);
        return -1;
    }

    int ret = 0;
    if (form == 3) {
        fuzz_write_fixed(fn, data, size);
    } else {
        static const char *mode[3] = {"wb", "wc", "wz"};
        samFile *out = sam_open(fn, mode[form]);
        if (!out) {
            ret = -1;
        } else {
            if (form == 1)
                hts_set_opt(out, CRAM_OPT_NO_REF, 1);
            bam1_t *b = bam_init1();
            if (!b)
                abort();
            if (sam_hdr_write(out, hdr) < 0)
                ret = -1;
            int r;
            while (ret == 0 && (r = sam_read1(in, hdr, b)) >= 0)
                if (sam_write1(out, hdr, b) < 0)
                    ret = -1;
            bam_destroy1(b);
            if (sam_close(out) < 0)
                ret = -1;
        }
    }
    hts_close(in);
    if (hdr_out)
        *hdr_out = hdr;
    else
        sam_hdr_destroy(hdr);
    return ret;
}

// Coordinate from one byte: near-zero, near-end, past-end or proportional
static hts_pos_t fuzz_pick_pos(uint8_t v, hts_pos_t len) {
    hts_pos_t p, k = v & 63;
    switch (v >> 6) {
    case 0:  p = k;            break;
    case 1:  p = len - k;      break;
    case 2:  p = len + k;      break;
    default: p = len * k / 63; break;
    }
    return p < 0 ? 0 : p;
}

static void fuzz_put_pos(kstring_t *ks, hts_pos_t p, int commas) {
    char buf[32];
    int n = snprintf(buf, sizeof(buf), "%"PRIhts_pos, p);
    for (int i = 0; i < n; i++) {
        if (commas && i > 0 && (n - i) % 3 == 0)
            kputc(',', ks);
        kputc(buf[i], ks);
    }
}

// Builds region string ks, and appends a matching BED line to bed
static void fuzz_make_region(kstring_t *ks, kstring_t *bed, const uint8_t *r,
                             sam_hdr_t *h) {
    int nref = h ? sam_hdr_nref(h) : 0;
    const char *name = "_none_";
    hts_pos_t len = 0;
    if (nref > 0) {
        int tid = ((r[0] >> 3) | ((r[3] >> 3) << 5)) % nref;
        name = sam_hdr_tid2name(h, tid);
        len = sam_hdr_tid2len(h, tid);
    }
    hts_pos_t beg = fuzz_pick_pos(r[1], len), end = fuzz_pick_pos(r[2], len);
    int commas = r[3] & 1, form = r[0] & 7;

    ks->l = 0;
    kputs("", ks); // ensure ks->s is allocated, even for an empty region
    switch (form) {
    case 0: kputs(name, ks); break;
    case 1: ksprintf(ks, "%s:", name); fuzz_put_pos(ks, beg, commas); break;
    case 2: case 3: case 6:
        ksprintf(ks, "%s:", name);
        fuzz_put_pos(ks, beg, commas);
        kputc('-', ks);
        if (form == 2) fuzz_put_pos(ks, end, commas);
        if (form == 6) kputs("99999999999999999999", ks);
        break;
    case 4:
        ksprintf(ks, "{%s}:", name);
        fuzz_put_pos(ks, beg, commas);
        kputc('-', ks);
        fuzz_put_pos(ks, end, commas);
        break;
    case 5: ksprintf(ks, "{%s}", name); break;
    default:
        switch ((r[3] >> 1) & 3) {
        case 0: kputs("*", ks); break;  // unplaced reads
        case 1: kputs(".", ks); break;  // whole file
        case 2:
            kputs("_absent_:", ks);
            fuzz_put_pos(ks, beg, commas);
            kputc('-', ks);
            fuzz_put_pos(ks, end, commas);
            break;
        case 3: break;                  // empty string
        }
        return; // no BED equivalent
    }

    // BED is 0-based, half open.  Whole-reference forms use 0..len.
    hts_pos_t bbeg = (form == 0 || form == 5) ? 0 : (beg > 0 ? beg - 1 : 0);
    hts_pos_t bend = (form == 0 || form == 5 || form == 3 || form == 1)
        ? len : end;
    ksprintf(bed, "%s\t%"PRIhts_pos"\t%"PRIhts_pos"\n", name, bbeg, bend);
}

static int fuzz_file_exists(const char *fn) {
    return access(fn, F_OK) == 0;
}

static void fuzz_file_reset(fuzz_file *f) {
    if (f->hdr) {
        sam_hdr_destroy(f->hdr);
        f->hdr = NULL;
    }
    free(f->extra_url);
    f->extra_url = NULL;
    f->extra_aln = 0;
}

static int fuzz_one_file_input(const uint8_t *data, size_t size,
                               size_t tool_bytes, int check_flags,
                               int flags,
                               const fuzz_step *steps, int nsteps) {
    fuzz_file *f = &fuzz_f;
    fuzz_begin_run(data, size);
    fuzz_file_reset(f); // anything left over from a run that exit()ed
    fuzz_extract_enter();
    fuzz_cleanup();

    if (setjmp(fuzz_jmp) != 0) {
        // exit() was called; don't add it to the test corpus.
        // Files and the header are tidied up at the start of the next run.
        return -1;
    }

    size_t ctrl = 1 + tool_bytes + FUZZ_MAX_REGIONS * FUZZ_REG_BYTES;
    if (size < ctrl + 5)
        return 0;

    const uint8_t *fb = data, *tool = data + 1, *reg = data + 1 + tool_bytes;
    data += ctrl;
    size -= ctrl;

    const uint8_t *extra = NULL;
    size_t extra_len = 0;
    if (flags & FUZZ_EXTRA) {
        if (size < 2)
            return 0;
        extra_len = data[0] | (data[1] << 8);
        if (extra_len > size - 2)
            extra_len = size - 2;
        extra = data + 2;
        data += 2 + extra_len;
        size -= 2 + extra_len;
        if (size < 5)
            return 0;
    }

    // Save the raw payloads first, in case the crash is in our preparation
    fuzz_extract_payload("fz.input", data, size);
    if (extra)
        fuzz_extract_payload("fz.extra", extra, extra_len);

    int is_aln = fuzz_check_input(data, size, check_flags & ~FUZZ_ANYFMT) == 0;
    if (!is_aln && !(check_flags & FUZZ_ANYFMT))
        return -1;

    static const char *ext[4] = {"bam", "cram", "sam.gz", "raw"};
    f->form = is_aln ? fb[0] & 3 : 3; // non-alignment data is always raw
    int idx_mode = (fb[0] >> 2) & 3;
    f->alt_idx = idx_mode == 2;
    snprintf(f->fn,       sizeof(f->fn),       "%s.%s", fuzz_prefix, ext[f->form]);
    snprintf(f->idx_fn,   sizeof(f->idx_fn),   "%s.xidx", fuzz_prefix);
    snprintf(f->idx_url,  sizeof(f->idx_url),  "%s##idx##%s", f->fn, f->idx_fn);
    snprintf(f->bed_fn,   sizeof(f->bed_fn),   "%s.bed", fuzz_prefix);
    snprintf(f->list_fn,  sizeof(f->list_fn),  "%s.list", fuzz_prefix);
    snprintf(f->extra_raw_fn, sizeof(f->extra_raw_fn), "%s.extra", fuzz_prefix);
    snprintf(f->extra_fn, sizeof(f->extra_fn), "%s.2.%s", fuzz_prefix, ext[f->form]);
    f->in = f->alt_idx ? f->idx_url : f->fn;

    if (!is_aln)
        fuzz_write_fixed(f->fn, data, size);   // f->hdr stays NULL
    else if (fuzz_write_aln(data, size, f->form, f->fn, &f->hdr) < 0)
        return -1;

    if (extra) {
        fuzz_write_fixed(f->extra_raw_fn, extra, extra_len);
        if (!(f->extra_url = fuzz_base64(extra, extra_len)))
            abort();
        if (extra_len >= 5 && fuzz_check_input(extra, extra_len, FUZZ_PARSE) == 0)
            f->extra_aln = fuzz_write_aln(extra, extra_len, f->form,
                                          f->extra_fn, NULL) == 0;
    }

    // BGZF SAM can only have a CSI index
    int min_shift = (idx_mode == 1 || f->form == 2) ? 14 : 0;
    if ((flags & FUZZ_INDEX) && idx_mode != 3) {
        (void) sam_index_build3(f->fn, f->alt_idx ? f->idx_fn : NULL,
                                min_shift, 0);
        if (f->extra_aln)
            (void) sam_index_build3(f->extra_fn, NULL, min_shift, 0);
    }

    // The index actually built, for -X (a missing one tests the error path)
    if (f->alt_idx) {
        snprintf(f->idx_path, sizeof(f->idx_path), "%s", f->idx_fn);
    } else {
        static const char *iext[3] = {"crai", "csi", "bai"};
        int i;
        for (i = 0; i < 3; i++) {
            snprintf(f->idx_path, sizeof(f->idx_path), "%s.%s", f->fn, iext[i]);
            if (fuzz_file_exists(f->idx_path))
                break;
        }
        if (i == 3)
            snprintf(f->idx_path, sizeof(f->idx_path), "%s.bai", f->fn);
    }

    // Regions, as strings and as a BED file
    kstring_t bed = KS_INITIALIZE;
    f->nreg = 1 + ((fb[0] >> 4) & 3);
    for (int i = 0; i < f->nreg; i++)
        fuzz_make_region(&f->reg[i], &bed, reg + i * FUZZ_REG_BYTES, f->hdr);
    fuzz_write_fixed(f->bed_fn, bed.s ? bed.s : "", bed.l);
    ks_free(&bed);

    // A file of filenames, for -b / -f style options
    char list[720];
    int ll = snprintf(list, sizeof(list), "%s\n", f->in);
    fuzz_write_fixed(f->list_fn, list, ll);

    // Build every step's command line first, so an extraction captures
    // them all even if the first one crashes
    static fuzz_args cmds[8];
    if (nsteps > 8)
        abort(); // harness bug
    for (int s = 0; s < nsteps; s++) {
        fuzz_args *a = &cmds[s];
        memset(a, 0, sizeof(*a));
        fuzz_arg(a, steps[s].main_fn ? "samtools" : steps[s].pseudo);
        steps[s].build(a, tool, f);
        fuzz_check_args(a, NULL);
    }
    if (*fuzz_extract_dir) {
        fuzz_write_repro(cmds, nsteps, NULL);
        if (fuzz_extract_only)
            return 0;
    }

    for (int s = 0; s < nsteps; s++) {
        if (!fuzz_quiet)
            fuzz_dump(cmds[s].argc, cmds[s].argv);
        if (!steps[s].main_fn) {
            fuzz_run_pseudo(&cmds[s]);
            continue;
        }

#ifdef __GLIBC__
        optind = 0;
#else
        optind = 1;
#endif
        opterr = 0;
        char **argv = fuzz_writable_argv(&cmds[s]);
        (void) steps[s].main_fn(cmds[s].argc - 1, argv + 1);
        fuzz_free_argv();
        fuzz_std_fds();
    }

    fuzz_file_reset(f);
    if (!*fuzz_extract_dir)
        fuzz_cleanup();
    return 0;
}

#endif /* FUZZ_COMMON_H */
