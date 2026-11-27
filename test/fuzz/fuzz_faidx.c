/*  test/fuzz/fuzz_faidx.c -- Fuzz driver for samtools faidx / fqidx

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

#include "fuzz_common.h"
#include <htslib/faidx.h>

int faidx_main(int argc, char *argv[]);
int fqidx_main(int argc, char *argv[]);

/*
Fuzzes samtools faidx and samtools fqidx

This source builds two fuzzers:

  Default:                 build-and-query
    The data is written to a file and indexed with "samtools faidx|fqidx",
    then queried with a second "samtools faidx|fqidx" call.  The index
    is always one samtools made itself, so this tests index construction
    and normal retrieval.

  -DFUZZ_SUPPLIED_INDEX:   query with a fuzzer-supplied index
    The fuzzer provides the .fai (and optionally .gzi) contents as well as
    the data, so we test loading of corrupt indices and of indices that
    don't match the data.  Kept as a separate binary because the input
    layout differs, so the corpora shouldn't be mixed.

Input layout:
    [OPT_BYTES option bytes][MAX_REGIONS * REG_BYTES region bytes]
    FUZZ_SUPPLIED_INDEX only: [fai_len 16-bit LE][gzi_len 16-bit LE]
                              [fai_len bytes .fai][gzi_len bytes .gzi]
    [sequence data: FASTA, FASTQ, or either BGZF compressed]

The region bytes pick sequence names and coordinates from the index (as
loaded by this harness), so most queries hit real sequences.  Coordinates
are chosen near 0, near the sequence end, past the end, or proportionally.

Build with ./build_faidx.sh (link faidx.o; add -DFUZZ_SUPPLIED_INDEX for
the second fuzzer).  See fuzz_common.h for general build and run details,
including FUZZ_EXTRACT for reproducing a crash outside the harness.  The
query regions are derived from the index the build step creates, so even
with FUZZ_EXTRACT_ONLY the build step runs; repro.sh is written before it
does, and rewritten with all the steps once the regions are known.

This fuzzer has to write files: the data, the indices, a regions file, and
optionally an -o output.  They use fuzz_prefix (_fz.<pid>.*), so parallel
jobs don't collide, and are removed at the start and end of every run.  Run
it from a scratch directory, ideally on tmpfs (e.g. /dev/shm), for speed.

Seed corpus files are made from small FASTA/FASTQ files with varied line
lengths.  19 below is CTRL_BYTES; all-zero control bytes is a valid choice.

Uncompressed, x.fa:
    samtools faidx x.fa                     # makes x.fa.fai
    (head -c 19 /dev/zero; cat x.fa) > seed_x                  # default
    (head -c 19 /dev/zero;
     perl -e 'print pack("vv", -s "x.fa.fai", 0)';
     cat x.fa.fai x.fa) > seed_x                               # supplied index

BGZF compressed.  Seed at least one of these: the fuzzer can't realistically
mutate plain text into valid BGZF (each block has a CRC), so without a seed
the .gzi code paths are never reached.
    bgzip -i -k x.fa                        # makes x.fa.gz and x.fa.gz.gzi
    samtools faidx x.fa.gz                  # makes x.fa.gz.fai
    (head -c 19 /dev/zero; cat x.fa.gz) > seed_xgz             # default
    (head -c 19 /dev/zero;
     perl -e 'print pack("vv", -s "x.fa.gz.fai", -s "x.fa.gz.gzi")';
     cat x.fa.gz.fai x.fa.gz.gzi x.fa.gz) > seed_xgz           # supplied index
The default fuzzer needs only the .gz, as samtools builds the .gzi itself.
The supplied-index fuzzer needs the .gzi too, as a compressed file can't be
queried without one.

No validity of result is assumed.  We're simply checking for crashes or
illegal program behaviour.
*/

#define OPT_BYTES   3
#define REG_BYTES   4
#define MAX_REGIONS 4
#define CTRL_BYTES  (OPT_BYTES + REG_BYTES * MAX_REGIONS)

/* ------------------------------------------------------------------------
 * Working files, all under fuzz_prefix
 */
static char fa_fn[320], fai_def[330], gzi_def[330], fai_alt[320],
    gzi_alt[330], reg_fn[320], out_fn[320], out_gz_fn[330];

static void init_filenames(void) {
    snprintf(fa_fn,     sizeof(fa_fn),     "%s.fa",      fuzz_prefix);
    snprintf(fai_def,   sizeof(fai_def),   "%s.fai",     fa_fn);
    snprintf(gzi_def,   sizeof(gzi_def),   "%s.gzi",     fa_fn);
    snprintf(fai_alt,   sizeof(fai_alt),   "%s.idx",     fuzz_prefix);
    snprintf(gzi_alt,   sizeof(gzi_alt),   "%s.idx.gzi", fuzz_prefix);
    snprintf(reg_fn,    sizeof(reg_fn),    "%s.reg",     fuzz_prefix);
    snprintf(out_fn,    sizeof(out_fn),    "%s.out",     fuzz_prefix);
    snprintf(out_gz_fn, sizeof(out_gz_fn), "%s.out.gz",  fuzz_prefix);
}

/* ------------------------------------------------------------------------
 * Region generation, using the shared coordinate pickers
 */

// r[0] bits 0-2: region syntax
// r[0] bits 3-7, r[3] bits 1-7: sequence selector
// r[1], r[2]: start and end coordinates
// r[3] bit 0: thousands separators in coordinates
static void make_region(kstring_t *ks, const uint8_t *r,
                        const char *name, hts_pos_t len) {
    hts_pos_t beg = fuzz_pick_pos(r[1], len);
    hts_pos_t end = fuzz_pick_pos(r[2], len);
    int commas = r[3] & 1;

    ks->l = 0;
    kputs("", ks); // ensure ks->s is allocated
    switch (r[0] & 7) {
    case 0: // name
        kputs(name, ks);
        break;
    case 1: // name:beg
        kputs(name, ks); kputc(':', ks); fuzz_put_pos(ks, beg, commas);
        break;
    case 2: // name:beg-end
        kputs(name, ks); kputc(':', ks); fuzz_put_pos(ks, beg, commas);
        kputc('-', ks); fuzz_put_pos(ks, end, commas);
        break;
    case 3: // name:beg-
        kputs(name, ks); kputc(':', ks); fuzz_put_pos(ks, beg, commas);
        kputc('-', ks);
        break;
    case 4: // {name}:beg-end, as used for names containing colons
        kputc('{', ks); kputs(name, ks); kputs("}:", ks);
        fuzz_put_pos(ks, beg, commas); kputc('-', ks);
        fuzz_put_pos(ks, end, commas);
        break;
    case 5: // {name}
        kputc('{', ks); kputs(name, ks); kputc('}', ks);
        break;
    case 6: // name:beg-<overflow>
        kputs(name, ks); kputc(':', ks); fuzz_put_pos(ks, beg, commas);
        kputs("-99999999999999999999", ks);
        break;
    default: // a name not in the index, for -c
        kputs("_absent_:", ks); fuzz_put_pos(ks, beg, commas);
        kputc('-', ks); fuzz_put_pos(ks, end, commas);
        break;
    }
}

/* ------------------------------------------------------------------------
 * Running samtools
 */

// Options which must be identical in the build and query steps
static void add_index_opts(fuzz_args *a, const uint8_t *opt) {
    if (opt[0] & 0x10)   ARG("-f");
    if (opt[2] & 0x01)   ARG2("--fai-idx", fai_alt);
    if (opt[2] & 0x02)   ARG2("--gzi-idx", gzi_alt);
}

// Runs a command: samtools, or the "rm" used between steps
static void run(fuzz_args *a, int (*tool_main)(int, char **)) {
    if (!fuzz_quiet)
        fuzz_dump(a->argc, a->argv);
    if (strcmp(a->argv[0], "rm") == 0) {
        for (int i = 2; i < a->argc; i++)  // skip "rm -f"
            unlink(a->argv[i]);
        return;
    }
#ifdef __GLIBC__
    optind = 0;
#else
    optind = 1;
#endif
    opterr = 0;
    char **argv = fuzz_writable_argv(a);
    (void) tool_main(a->argc - 1, argv + 1);
    fuzz_free_argv();
}

/*
 * Option bytes:
 *
 * opt[0] bits 0-1  number of regions - 1
 *        bit  2    pass regions via -r FILE instead of on the command line
 *        bit  3    use fqidx instead of faidx
 *        bit  4    -f
 *        bit  5    unused
 *        bit  6    delete the index before querying, so the query step
 *                  rebuilds it (build-and-query fuzzer only)
 *        bit  7    -c
 * opt[1] bit  0    -i
 *        bits 1-2  --mark-strand
 *        bits 3-4  -n
 *        bit  5    -o FILE
 *        bit  6    -o FILE is .gz (BGZF output)
 *        bit  7    --output-fmt-option level=1
 * opt[2] bit  0    --fai-idx with a non-default name
 *        bit  1    --gzi-idx with a non-default name
 *        bits 2-7  unused
 */
int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    static kstring_t region[MAX_REGIONS];
    static fuzz_args cmds[3];   // build, rm, query
    int ncmds = 0;

    fuzz_extract_enter();
    init_filenames();
    fuzz_cleanup();

    if (setjmp(fuzz_jmp) != 0) {
        // exit() was called; don't add it to the test corpus.
        // Files are tidied up at the start of the next run.
        return -1;
    }

    if (size < CTRL_BYTES)
        return 0;

    const uint8_t *opt = data, *reg = data + OPT_BYTES;
    data += CTRL_BYTES;
    size -= CTRL_BYTES;
    fuzz_extract_payload("fz.input", data, size);

    int fastq_tool = opt[0] & 0x08;
    int (*tool_main)(int, char **) = fastq_tool ? fqidx_main : faidx_main;
    const char *tool   = fastq_tool ? "fqidx" : "faidx";
    const char *fai_fn = (opt[2] & 0x01) ? fai_alt : fai_def;
    const char *gzi_fn = (opt[2] & 0x02) ? gzi_alt : gzi_def;
    fuzz_args *a;
    int already_run = 0; // commands in cmds[] that have been run already

#ifdef FUZZ_SUPPLIED_INDEX
    // Step 1: write the fuzzer-supplied index files and the data
    if (size < 4)
        return 0;
    size_t fai_len = data[0] | (data[1] << 8);
    size_t gzi_len = data[2] | (data[3] << 8);
    data += 4;
    size -= 4;

    if (fai_len > size)
        fai_len = size;
    fuzz_write_fixed(fai_fn, data, fai_len);
    data += fai_len;
    size -= fai_len;

    if (gzi_len > size)
        gzi_len = size;
    if (gzi_len)
        fuzz_write_fixed(gzi_fn, data, gzi_len);
    data += gzi_len;
    size -= gzi_len;

    fuzz_write_fixed(fa_fn, data, size);
#else
    // Step 1: write the data and index it
    fuzz_write_fixed(fa_fn, data, size);

    a = memset(&cmds[ncmds++], 0, sizeof(*cmds));
    ARG("samtools");
    ARG(tool);
    add_index_opts(a, opt);
    ENDOPTS();
    ARG(fa_fn);
    fuzz_check_args(a, NULL);
    if (*fuzz_extract_dir)
        fuzz_write_repro(cmds, ncmds, NULL);
    run(a, tool_main);
    already_run = 1;
#endif

    // Step 2: choose regions using the sequence names and lengths in the
    // index.  If it won't load, the regions use a name that doesn't exist.
    enum fai_format_options fmt =
        (fastq_tool || (opt[0] & 0x10)) ? FAI_FASTQ : FAI_FASTA;
    faidx_t *fai = fai_load3_format(fa_fn, fai_fn, gzi_fn, 0, fmt);
    if (!fai)
        fai = fai_load3_format(fa_fn, fai_fn, gzi_fn, 0,
                               fmt == FAI_FASTQ ? FAI_FASTA : FAI_FASTQ);
    int nseq = fai ? faidx_nseq(fai) : 0;
    int nreg = 1 + (opt[0] & 3);

    for (int i = 0; i < nreg; i++) {
        const uint8_t *r = reg + i * REG_BYTES;
        const char *name = "_none_";
        hts_pos_t len = 0;
        if (nseq > 0) {
            int idx = ((r[0] >> 3) | ((r[3] >> 1) << 5)) % nseq;
            const char *n = faidx_iseq(fai, idx);
            if (n) {
                name = n;
                len = faidx_seq_len(fai, name);
                if (len < 0)
                    len = 0;
            }
        }
        make_region(&region[i], r, name, len);
    }
    if (fai)
        fai_destroy(fai);

#ifndef FUZZ_SUPPLIED_INDEX
    // Optionally force the query step to rebuild the index itself
    if (opt[0] & 0x40) {
        a = memset(&cmds[ncmds++], 0, sizeof(*cmds));
        ARG("rm");
        ARG("-f");
        ARG(fai_fn);
        ARG(gzi_fn);
    }
#endif

    // Step 3: query
    if (opt[0] & 0x04) {
        kstring_t rf = KS_INITIALIZE;
        for (int i = 0; i < nreg; i++) {
            kputs(region[i].s, &rf);
            kputc('\n', &rf);
        }
        fuzz_write_fixed(reg_fn, rf.s, rf.l);
        ks_free(&rf);
    }

    a = memset(&cmds[ncmds++], 0, sizeof(*cmds));
    ARG("samtools");
    ARG(tool);
    add_index_opts(a, opt);
    if (opt[0] & 0x80)   ARG("-c");
    if (opt[1] & 0x01)   ARG("-i");
    CHOOSE("--mark-strand", (opt[1] >> 1) & 3, NULL, "no", "sign", "custom, fwd,_rev");
    CHOOSE("-n",            (opt[1] >> 3) & 3, NULL, "0",  "1",    "7");
    if (opt[1] & 0x20)   ARG2("-o", (opt[1] & 0x40) ? out_gz_fn : out_fn);
    if (opt[1] & 0x80)   ARG2("--output-fmt-option", "level=1");
    // Threads make runs non-deterministic.
    //if (opt[2] & 0x04) ARG2("-@", "2");
    if (opt[0] & 0x04)   ARG2("-r", reg_fn);

    // Sequence names, and hence regions, may start with '-'.  GNU getopt
    // permutes arguments, so without "--" a region such as "-ofoo:1-10"
    // would be parsed as "-o foo:1-10" and write a file.
    ENDOPTS();
    ARG(fa_fn);
    if (!(opt[0] & 0x04))
        for (int i = 0; i < nreg; i++)
            ARG_DATA(region[i].s);
    fuzz_check_args(a, NULL);

    if (*fuzz_extract_dir) {
        fuzz_write_repro(cmds, ncmds, NULL);
        if (fuzz_extract_only)
            return 0;
    }
    for (int c = already_run; c < ncmds; c++)
        run(&cmds[c], tool_main);

    if (!*fuzz_extract_dir)
        fuzz_cleanup();
    return 0;
}
