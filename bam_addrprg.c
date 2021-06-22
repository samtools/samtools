/* bam_addrprg.c -- samtools command to add or replace readgroups.

   Copyright (c) 2013, 2015-2017, 2019-2021 Genome Research Limited.

   Author: Martin O. Pollard <mp15@sanger.ac.uk>

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

#include <config.h>

#include <htslib/sam.h>
#include <htslib/kstring.h>
#include "samtools.h"
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

typedef enum {
    overwrite_all,
    orphan_only,
} rg_mode;

struct parsed_opts {
    char* input_name;
    char* output_name;
    char* rg_id;
    char* rg_line;
    int no_pg;
    rg_mode mode;
    sam_global_args ga;
    htsThreadPool p;
    int uncompressed;
    int overwrite_hdr_rg;
};

struct state;
typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

struct state {
    samFile* input_file;
    sam_hdr_t* input_header;
    samFile* output_file;
    sam_hdr_t* output_header;
    char* rg_id;
    void (*mode_func)(const state_t*, bam1_t*);
};

static void cleanup_opts(parsed_opts_t* opts)
{
    if (!opts) return;
    free(opts->rg_id);
    free(opts->output_name);
    free(opts->input_name);
    free(opts->rg_line);
    if (opts->p.pool) hts_tpool_destroy(opts->p.pool);
    sam_global_args_free(&opts->ga);
    free(opts);
}

static void cleanup_state(state_t* state)
{
    if (!state) return;
    free(state->rg_id);
    if (state->output_file) sam_close(state->output_file);
    sam_hdr_destroy(state->output_header);
    if (state->input_file) sam_close(state->input_file);
    sam_hdr_destroy(state->input_header);
    free(state);
}

// Converts \t and \n into real tabs and newlines
static char* basic_unescape(const char* in)
{
    assert(in);
    char *ptr, *out;
    out = ptr = malloc(strlen(in)+1);
    size_t size = 0;
    while (*in) {
        if (*in == '\\') {
            ++in;
            if (*in == '\0') {
                fprintf(stderr, "[%s] Unterminated escape sequence.\n", __func__);
                free(out);
                return NULL;
            }
            switch (*in) {
            case '\\':
                *ptr = '\\';
                break;
            case 't':
                *ptr = '\t';
                break;
            case 'n':
                fprintf(stderr, "[%s] \\n in escape sequence is not supported.\n", __func__);
                free(out);
                return NULL;
            default:
                fprintf(stderr, "[%s] Unsupported escape sequence.\n", __func__);
                free(out);
                return NULL;
            }
        } else {
            *ptr = *in;
        }
        ++in;
        ++ptr;
        ++size;
    }
    *ptr = '\0';
    ++size;
    char* tmp = (char*)realloc(out, size);
    if (!tmp) {
        free(out);
    }
    return tmp;
}

// Malloc a string containing [s,slim) or to the end of s if slim is NULL.
// If lenp is non-NULL, stores the length of the resulting string there.
static char *dup_substring(const char *s, const char *slim, size_t *lenp)
{
    size_t len = slim? (slim - s) : strlen(s);
    char *ns = malloc(len+1);
    if (ns == NULL) return NULL;
    memcpy(ns, s, len);
    ns[len] = '\0';
    if (lenp) *lenp = len;
    return ns;
}


// Given a @RG line return the id
static char* get_rg_id(const char *line)
{
    const char *id = strstr(line, "\tID:");
    if (! id) return NULL;

    id += 4;
    return dup_substring(id, strchr(id, '\t'), NULL);
}


static void usage(FILE *fp)
{
    fprintf(fp,
            "Usage: samtools addreplacerg [options] [-r <@RG line> | -R <existing id>] [-m orphan_only|overwrite_all] [-o <output.bam>] <input.bam>\n"
            "\n"
            "Options:\n"
            "  -m MODE   Set the mode of operation from one of overwrite_all, orphan_only [overwrite_all]\n"
            "  -o FILE   Where to write output to [stdout]\n"
            "  -r STRING @RG line text\n"
            "  -R STRING ID of @RG line in existing header to use\n"
            "  -u        Output uncompressed data\n"
            "  -w        Overwrite an existing @RG line\n"
            "  --no-PG   Do not add a PG line\n"
            );
    sam_global_opt_help(fp, "..O..@..");
}

static bool parse_args(int argc, char** argv, parsed_opts_t** opts)
{
    *opts = NULL;
    int n;

    if (argc == 1) { usage(stdout); return true; }

    parsed_opts_t* retval = calloc(1, sizeof(parsed_opts_t));
    if (! retval ) {
        fprintf(stderr, "[%s] Out of memory allocating parsed_opts_t\n", __func__);
        return false;
    }
    // Set defaults
    retval->mode = overwrite_all;
    sam_global_args_init(&retval->ga);
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS(0, 0, 'O', 0, 0, '@'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };
    kstring_t rg_line = {0,0,NULL};

    while ((n = getopt_long(argc, argv, "r:R:m:o:O:h@:uw", lopts, NULL)) >= 0) {
        switch (n) {
            case 'r':
                // Are we adding to existing rg line?
                if (ks_len(&rg_line) == 0) {
                    if (strlen(optarg)<3 || (optarg[0] != '@' && optarg[1] != 'R' && optarg[2] != 'G')) {
                        kputs("@RG\t", &rg_line);
                    }
                } else {
                    kputs("\t", &rg_line);
                }
                kputs(optarg, &rg_line);
                break;
            case 'R':
                retval->rg_id = strdup(optarg);
                break;
            case 'm': {
                if (strcmp(optarg, "overwrite_all") == 0) {
                    retval->mode = overwrite_all;
                } else if (strcmp(optarg, "orphan_only") == 0) {
                    retval->mode = orphan_only;
                } else {
                    usage(stderr);
                    return false;
                }
                break;
            }
            case 'o':
                retval->output_name = strdup(optarg);
                break;
            case 'h':
                usage(stdout);
                free(retval);
                return true;
            case 1:
                retval->no_pg = 1;
                break;
            case 'u':
                retval->uncompressed = 1;
                break;
            case 'w':
                retval->overwrite_hdr_rg = 1;
                break;
            case '?':
                usage(stderr);
                free(retval);
                return false;
            case 'O':
            default:
                if (parse_sam_global_opt(n, optarg, lopts, &retval->ga) == 0) break;
                usage(stderr);
                free(retval);
                return false;
        }
    }
    retval->rg_line = ks_release(&rg_line);

    if (argc-optind < 1) {
        fprintf(stderr, "You must specify an input file.\n");
        usage(stderr);
        cleanup_opts(retval);
        return false;
    }
    if (retval->rg_id && retval->rg_line) {
        fprintf(stderr, "The options -r and -R are mutually exclusive.\n");
        cleanup_opts(retval);
        return false;
    }

    if (retval->rg_line)
    {
        char* tmp = basic_unescape(retval->rg_line);

        if ((retval->rg_id = get_rg_id(tmp)) == NULL) {
            fprintf(stderr, "[%s] The supplied RG line lacks an ID tag.\n", __func__);
            free(tmp);
            cleanup_opts(retval);
            return false;
        }
        free(retval->rg_line);
        retval->rg_line = tmp;
    }
    retval->input_name = strdup(argv[optind+0]);

    if (retval->ga.nthreads > 0) {
        if (!(retval->p.pool = hts_tpool_init(retval->ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            return false;
        }
    }

    *opts = retval;
    return true;
}

static void overwrite_all_func(const state_t* state, bam1_t* file_read)
{
    uint8_t* data = (uint8_t*)strdup(state->rg_id);
    int len = strlen(state->rg_id)+1;
    // If the old exists delete it
    uint8_t* old = bam_aux_get(file_read, "RG");
    if (old != NULL) {
        bam_aux_del(file_read, old);
    }

    bam_aux_append(file_read, "RG", 'Z', len, data);
    free(data);
}

static void orphan_only_func(const state_t* state, bam1_t* file_read)
{
    uint8_t* data = (uint8_t*)strdup(state->rg_id);
    int len = strlen(state->rg_id)+1;
    // If the old exists don't do anything
    uint8_t* old = bam_aux_get(file_read, "RG");
    if (old == NULL) {
        bam_aux_append(file_read, "RG",'Z',len,data);
    }
    free(data);
}

static bool init(const parsed_opts_t* opts, state_t** state_out) {
    char output_mode[9] = "w";
    state_t* retval = (state_t*) calloc(1, sizeof(state_t));

    if (retval == NULL) {
        fprintf(stderr, "[init] Out of memory allocating state struct.\n");
        return false;
    }
    *state_out = retval;

    // Open files
    retval->input_file = sam_open_format(opts->input_name, "r", &opts->ga.in);
    if (retval->input_file == NULL) {
        print_error_errno("addreplacerg", "could not open \"%s\"", opts->input_name);
        return false;
    }
    retval->input_header = sam_hdr_read(retval->input_file);

    retval->output_header = sam_hdr_dup(retval->input_header);

    if (opts->uncompressed)
        strcat(output_mode, "0");
    if (opts->output_name) // File format auto-detection
        sam_open_mode(output_mode + strlen(output_mode),
                      opts->output_name, NULL);
    retval->output_file = sam_open_format(opts->output_name == NULL?"-":opts->output_name, output_mode, &opts->ga.out);

    if (retval->output_file == NULL) {
        print_error_errno("addreplacerg", "could not create \"%s\"", opts->output_name);
        return false;
    }

    if (opts->p.pool) {
        hts_set_opt(retval->input_file,  HTS_OPT_THREAD_POOL, &opts->p);
        hts_set_opt(retval->output_file, HTS_OPT_THREAD_POOL, &opts->p);
    }

    if (opts->rg_line) {
        // Append new RG line to header.
        // Check does not already exist
        kstring_t hdr_line = { 0, 0, NULL };
        if (sam_hdr_find_line_id(retval->output_header, "RG", "ID", opts->rg_id, &hdr_line) == 0) {
            if (opts->overwrite_hdr_rg) {
                if(-1 == sam_hdr_remove_line_id(retval->output_header, "RG", "ID", opts->rg_id)) {
                    fprintf(stderr, "[init] Error removing the RG line with ID:%s from the output header.\n", opts->rg_id);
                    ks_free(&hdr_line);
                    return false;
                }
            } else {
                fprintf(stderr, "[init] RG line with ID:%s already present in the header. Use -w to overwrite.\n", opts->rg_id);
                ks_free(&hdr_line);
                return false;
            }
        }
        ks_free(&hdr_line);

        if (-1 == sam_hdr_add_lines(retval->output_header, opts->rg_line, strlen(opts->rg_line))) {
            fprintf(stderr, "[init] Error adding RG line with ID:%s to the output header.\n", opts->rg_id);
            return false;
        }
        if (opts->mode == overwrite_all &&
            -1 == sam_hdr_remove_except(retval->output_header, "RG", "ID", opts->rg_id)) {
            fprintf(stderr, "[init] Error removing the old RG lines from the output header.\n");
            return false;
        }
        retval->rg_id = strdup(opts->rg_id);
    } else {
        if (opts->rg_id) {
            // Confirm what has been supplied exists
            kstring_t hdr_line = { 0, 0, NULL };
            if (sam_hdr_find_line_id(retval->output_header, "RG", "ID", opts->rg_id, &hdr_line) < 0) {
                fprintf(stderr, "RG ID supplied does not exist in header. Supply full @RG line with -r instead?\n");
                return false;
            }
            retval->rg_id = strdup(opts->rg_id);
            ks_free(&hdr_line);
        } else {
            kstring_t rg_id = { 0, 0, NULL };
            if (sam_hdr_find_tag_id(retval->output_header, "RG", NULL, NULL, "ID", &rg_id) < 0) {
                fprintf(stderr, "No RG specified on command line or in existing header.\n");
                return false;
            }
            retval->rg_id = ks_release(&rg_id);
        }
    }

    switch (opts->mode) {
        case overwrite_all:
            retval->mode_func = &overwrite_all_func;
            break;
        case orphan_only:
            retval->mode_func = &orphan_only_func;
            break;
    }

    return true;
}

static bool readgroupise(parsed_opts_t *opts, state_t* state, char *arg_list)
{
    if (!opts->no_pg && sam_hdr_add_pg(state->output_header, "samtools",
                                       "VN", samtools_version(),
                                       arg_list ? "CL": NULL,
                                       arg_list ? arg_list : NULL,
                                       NULL))
        return false;

    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        print_error_errno("addreplacerg", "[%s] Could not write header to output file", __func__);
        return false;
    }
    char *idx_fn = NULL;
    if (opts->ga.write_index) {
        if (!(idx_fn = auto_index(state->output_file, opts->output_name, state->output_header)))
            return false;
    }

    bam1_t* file_read = bam_init1();
    int ret;
    while ((ret = sam_read1(state->input_file, state->input_header, file_read)) >= 0) {
        state->mode_func(state, file_read);

        if (sam_write1(state->output_file, state->output_header, file_read) < 0) {
            print_error_errno("addreplacerg", "[%s] Could not write read to output file", __func__);
            bam_destroy1(file_read);
            free(idx_fn);
            return false;
        }
    }
    bam_destroy1(file_read);
    if (ret != -1) {
        print_error_errno("addreplacerg", "[%s] Error reading from input file", __func__);
        free(idx_fn);
        return false;
    } else {

        if (opts->ga.write_index) {
            if (sam_idx_save(state->output_file) < 0) {
                print_error_errno("addreplacerg", "[%s] Writing index failed", __func__);
                free(idx_fn);
                return false;
            }
        }
        free(idx_fn);
        return true;
    }
}

int main_addreplacerg(int argc, char** argv)
{
    parsed_opts_t* opts = NULL;
    state_t* state = NULL;
    char *arg_list = stringify_argv(argc+1, argv-1);
    if (!arg_list)
        return EXIT_FAILURE;

    if (!parse_args(argc, argv, &opts)) goto error;
    if (opts) { // Not an error but user doesn't want us to proceed
        if (!init(opts, &state) || !readgroupise(opts, state, arg_list))
            goto error;
    }

    cleanup_state(state);
    cleanup_opts(opts);
    free(arg_list);

    return EXIT_SUCCESS;
error:
    cleanup_state(state);
    cleanup_opts(opts);
    free(arg_list);

    return EXIT_FAILURE;
}
