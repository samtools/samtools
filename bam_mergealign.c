/* bam_mergealign.c -- samtools command to merge alignments into another CRAM.

   Copyright (c) 2024 Genome Research Limited.

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
#include <htslib/thread_pool.h>
#include "samtools.h"
#include "sam_opts.h"


#include <stdbool.h>

typedef enum {
    baseline,
    minimap2,
    ngmlr,
    bwa_mem
} aligner_profile;

struct parsed_opts {
    char* arg_list;
    char* input_base_name;
    char* input_align_name;
    char* output_name;
    aligner_profile profile;
    sam_global_args ga;
    htsThreadPool p;
    int uncompressed;
    int no_pg;
};

struct state;
typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

struct state {
    samFile* input_base_file;
    sam_hdr_t* input_base_header;
    samFile* input_align_file;
    sam_hdr_t* input_align_header;
    samFile* output_file;
    sam_hdr_t* output_header;
    char* rg_id;
    void (*mode_func)(const state_t*, bam1_t*);
};

/////////////////////////////////////////////////

static void cleanup_opts(parsed_opts_t* opts)
{
    if (!opts) return;
    free(opts->output_name);
    free(opts->input_align_name);
    free(opts->input_base_name);
    if (opts->p.pool) hts_tpool_destroy(opts->p.pool);
    sam_global_args_free(&opts->ga);
    free(opts);
}

static void cleanup_state(state_t* state)
{
    if (!state) return;

    if (state->output_file) sam_close(state->output_file);
    sam_hdr_destroy(state->output_header);
    if (state->input_align_file) sam_close(state->input_align_file);
    if (state->input_base_file) sam_close(state->input_base_file);
    sam_hdr_destroy(state->input_align_header);
    sam_hdr_destroy(state->input_base_header);
    free(state);
}

/////////////////////////////////////////////

static bool init(const parsed_opts_t* opts, state_t** state_out) {
    char output_mode[9] = "w";
    state_t* retval = (state_t*) calloc(1, sizeof(state_t));

    if (retval == NULL) {
        fprintf(stderr, "[merge_align:init] Out of memory allocating state struct.\n");
        return false;
    }
    *state_out = retval;

    // Open files
    // align
    retval->input_align_file = sam_open_format(opts->input_align_name, "r", &opts->ga.in);
    if (retval->input_align_file == NULL) {
        print_error_errno("merge_align", "could not open \"%s\"", opts->input_align_name);
        return false;
    }
    retval->input_align_header = sam_hdr_read(retval->input_align_file);
    // base
    retval->input_base_file = sam_open_format(opts->input_base_name, "r", &opts->ga.in);
    if (retval->input_base_file == NULL) {
        print_error_errno("merge_align", "could not open \"%s\"", opts->input_base_name);
        return false;
    }
    retval->input_align_header = sam_hdr_read(retval->input_align_file);

    retval->output_header = sam_hdr_dup(retval->input_align_header);

    if (opts->uncompressed)
        strncat(output_mode, "0", 9);
    if (opts->output_name) // File format auto-detection
        sam_open_mode(output_mode + strlen(output_mode),
                      opts->output_name, NULL);
    retval->output_file = sam_open_format(opts->output_name == NULL?"-":opts->output_name, output_mode, &opts->ga.out);

    if (retval->output_file == NULL) {
        print_error_errno("merge_align", "could not create \"%s\"", opts->output_name);
        return false;
    }

    if (opts->p.pool) {
        hts_set_opt(retval->input_base_file,  HTS_OPT_THREAD_POOL, &opts->p);
        hts_set_opt(retval->input_align_file,  HTS_OPT_THREAD_POOL, &opts->p);
        hts_set_opt(retval->output_file, HTS_OPT_THREAD_POOL, &opts->p);
    }
    return true;
}

/////////////////////////////////////////////

static void usage(FILE *fp)
{
    fprintf(fp,
            "Usage: samtools mergealign [options] [-o <output.bam>] <input_base.bam> <input_aligner.bam>\n"
            "\n"
            "Options:\n"
            "  -p PROFILE Aligner profile to apply to tag copying [baseline]\n"
            "  -o FILE    Where to write output to [stdout]\n"
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

    retval->arg_list = stringify_argv(argc+1, argv-1);
    if (!retval->arg_list)
        return false;

    // Set defaults
    retval->profile = baseline;
    sam_global_args_init(&retval->ga);
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS(0, 0, 'O', 0, 0, '@'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    while ((n = getopt_long(argc, argv, "p:o:O:h@:u", lopts, NULL)) >= 0) {
        switch (n) {
            case 'p': {
                if (strcmp(optarg, "baseline") == 0) {
                    retval->profile = baseline;
                } else if (strcmp(optarg, "minimap2") == 0) {
                    retval->profile = minimap2;
                } else if (strcmp(optarg, "ngmlr") == 0) {
                    retval->profile = ngmlr;
                } if (strcmp(optarg, "bwa_mem") == 0) {
                    retval->profile = bwa_mem;
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

    if (argc-optind < 2) {
        fprintf(stderr, "You must specify two input files.\n");
        usage(stderr);
        cleanup_opts(retval);
        return false;
    }
    retval->input_base_name = strdup(argv[optind+0]);
    retval->input_align_name = strdup(argv[optind+1]);

    if (retval->ga.nthreads > 0) {
        if (!(retval->p.pool = hts_tpool_init(retval->ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            return false;
        }
    }

    *opts = retval;
    return true;
}

/////////////////////////////////////////////////

static void bam_readid(char* const name_buff, const bam1_t* read)
{
    snprintf(name_buff, 255, "%.254s%c", bam_get_qname(read), ((read)->core.flag & (BAM_FREAD1|BAM_FREAD2) << 6)+'0');
}

// Assumptions: inputs are consistently sorted
// Read a template read from input_base (skip all non primary)
// POINT B
// Read a working read from input_align
// key = readname + readnumber
// while working_key != base_key read a template read from input_base_file
// POINT A
// Strip banned tags (SA, MD, MM) from working read
// Copy all non-special AUX to working read
// Clip MM tag and copy to working read
// write working read to output
// GOTO B
static bool mergealign(const parsed_opts_t *opts, state_t* state)
{
    if (!opts->no_pg && sam_hdr_add_pg(state->output_header, "samtools",
                                       "VN", samtools_version(),
                                       opts->arg_list ? "CL": NULL,
                                       opts->arg_list ? opts->arg_list : NULL,
                                       NULL))
        return false;

    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        print_error_errno("mergealign", "[%s] Could not write header to output file", __func__);
        return false;
    }
    char *idx_fn = NULL;
    if (opts->ga.write_index) {
        if (!(idx_fn = auto_index(state->output_file, opts->output_name, state->output_header)))
            return false;
    }

    bam1_t* file_base_read = bam_init1();
    bam1_t* file_align_read = bam_init1();
    int ret;
    char align_name_buff[255];
    char base_name_buff[255];
    if ((ret = sam_read1(state->input_base_file, state->input_base_header, file_base_read)) < 0) {
        print_error_errno("mergealign", "[%s] Error reading from input file", __func__);
        free(idx_fn);
        return false;
    }
    bam_readid(base_name_buff, file_base_read);

    while ((ret = sam_read1(state->input_align_file, state->input_align_header, file_align_read)) >= 0) {
        bam_readid(align_name_buff, file_align_read);
        // Read from base file until we get a matching readname + READ1 READ2 flag
        while (!strcmp(base_name_buff, align_name_buff)) {
            if ((ret = sam_read1(state->input_base_file, state->input_base_header, file_base_read)) < 0) {
                print_error_errno("mergealign", "[%s] Error reading from input file", __func__);
                bam_destroy1(file_align_read);
                bam_destroy1(file_base_read);
                free(idx_fn);
                return false;
            }
            // This shouldn't happen if both are name sorted, ba
            if (ret == 0) {
                print_error_errno("mergealign", "[%s] Error, inputs from base file exhausted", __func__);
                bam_destroy1(file_align_read);
                bam_destroy1(file_base_read);
                free(idx_fn);
                return false;
            }
            // Skip secondary and supplementary reads in basefile
            if (file_base_read->core.flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY) != 0) {
                continue;
            }
            bam_readid(base_name_buff, file_base_read);
        }
        // Remove banned tags
        for (uint8_t* aux = bam_aux_first(file_align_read); aux != NULL; ) {
            if (!strncmp(bam_aux_tag(aux), "SA",2)) {
                aux = bam_aux_remove(file_align_read, aux);
            } else {
                aux = bam_aux_next(file_align_read, aux);
            }
        }
        // Copy specified tags
        const char copy_list[][2] = {"SA", "MD"};
        uint8_t* copy_tag;
        for (int i = 0; i < sizeof(copy_list); i++){
            const char* copy_tag_id = copy_list[i];
            // Does source read have the tag?
            if(( copy_tag = bam_aux_get(file_base_read, copy_tag_id)) != NULL) {
                // Yes copy it
                if (bam_aux_append_from_aux(file_align_read, file_base_read, copy_tag) == -1) {
                    print_error_errno("mergealign", "[%s] Error copying aux tags to record", __func__);
                    bam_destroy1(file_align_read);
                    bam_destroy1(file_base_read);
                    free(idx_fn);
                    return false;
                }
            }
        }
        // check nothing went horribly wrong
        if (errno != ENOENT) {
                print_error_errno("mergealign", "[%s] Error removing aux tags from record", __func__);
                bam_destroy1(file_align_read);
                bam_destroy1(file_base_read);
                free(idx_fn);
                return false;
        }
        // Write out result
        if (sam_write1(state->output_file, state->output_header, file_align_read) < 0) {
            print_error_errno("mergealign", "[%s] Could not write read to output file", __func__);
            bam_destroy1(file_align_read);
            bam_destroy1(file_base_read);
            free(idx_fn);
            return false;
        }
    }
    bam_destroy1(file_align_read);
    bam_destroy1(file_base_read);
    if (ret != -1) {
        print_error_errno("mergealign", "[%s] Error reading from input file", __func__);
        free(idx_fn);
        return false;
    } else {
        if (opts->ga.write_index) {
            if (sam_idx_save(state->output_file) < 0) {
                print_error_errno("mergealign", "[%s] Writing index failed", __func__);
                free(idx_fn);
                return false;
            }
        }
        free(idx_fn);
        return true;
    }
}

int main_merge_align(int argc, char ** argv)
{
    parsed_opts_t* opts = NULL;
    state_t* state = NULL;
    if (!parse_args(argc, argv, &opts)) goto error;
    if (opts) { // Not an error but user doesn't want us to proceed
        if (!init(opts, &state) || !mergealign(opts, state))
            goto error;
    }

    cleanup_state(state);
    cleanup_opts(opts);
    return EXIT_SUCCESS;
error:
    cleanup_state(state);
    cleanup_opts(opts);

    return EXIT_FAILURE;
}