/*  bam_split.c -- split subcommand.

    Copyright (C) 2013-2016,2018-2019,2023,2024 Genome Research Ltd.

    Author: Martin Pollard <mp15@sanger.ac.uk>

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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <unistd.h>
#include <regex.h>
#include <assert.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>
#include <htslib/cram.h>
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include "samtools.h"


KHASH_MAP_INIT_STR(c2i, int)

struct parsed_opts {
    const char *merged_input_name;
    const char *unaccounted_header_name;
    const char *unaccounted_name;
    const char *output_format_string;
    const char *tag;
    long max_split;
    bool verbose;
    int no_pg;
    int zero_pad;
    sam_global_args ga;
};

#define DEFAULT_MAX_SPLIT 100

typedef struct parsed_opts parsed_opts_t;

struct state {
    samFile* merged_input_file;
    sam_hdr_t* merged_input_header;
    samFile* unaccounted_file;
    sam_hdr_t* unaccounted_header;
    char *unaccounted_idx_fn;
    size_t output_count;
    char **tag_vals;
    char **index_file_name;
    char **output_file_name;
    samFile **output_file;
    sam_hdr_t **output_header;
    kh_c2i_t* tag_val_hash;
    htsThreadPool p;
    int write_index;
};

typedef struct state state_t;

static int cleanup_state(state_t* status, bool check_close);
static void cleanup_opts(parsed_opts_t* opts);

static void usage(FILE *write_to)
{
    fprintf(write_to,
"Usage: samtools split [-u <unaccounted.bam>] [-h <unaccounted_header.sam>]\n"
"                      [-f <format_string>] [-v] <merged.bam>\n"
"Options:\n"
"  -f STRING           output filename format string [\"%%*_%%#.%%.\"]\n"
"  -u FILE1            put left-over reads in FILE1\n"
"  -h FILE2            ... and override the header with FILE2 (-u file only)\n"
"  -d TAG              split by TAG value. TAG value must be a string.\n"
"  -p NUMBER           zero-pad numbers in filenames to NUMBER digits\n"
"  -M,--max-split NUM  limit number of output files from -d to NUM [%d]\n"
"  -v                  verbose output\n"
"  --no-PG             do not add a PG line\n",
            DEFAULT_MAX_SPLIT);
    sam_global_opt_help(write_to, "-....@..");
    fprintf(write_to,
"\n"
"Format string expansions:\n"
"  %%%%     %%\n"
"  %%*     basename\n"
"  %%#     index (of @RG in the header, or count of TAG values seen so far)\n"
"  %%!     @RG ID or TAG value\n"
"  %%.     filename extension for output format\n"
      );
}

// Takes the command line options and turns them into something we can understand
static parsed_opts_t* parse_args(int argc, char** argv)
{
    if (argc == 1) { usage(stdout); return NULL; }

    const char *optstring = "vf:h:u:d:M:p:@:";
    char *default_format_string = "%*_%#.%.";

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0, '@'),
        {"no-PG", no_argument, NULL, 1},
        {"max-split", required_argument, NULL, 'M'},
        {"zero-pad", required_argument, NULL, 'p'},
        { NULL, 0, NULL, 0 }
    };

    parsed_opts_t* retval = calloc(sizeof(parsed_opts_t), 1);
    if (! retval ) { perror("cannot allocate option parsing memory"); return NULL; }

    retval->max_split = DEFAULT_MAX_SPLIT;
    sam_global_args_init(&retval->ga);

    int opt;
    while ((opt = getopt_long(argc, argv, optstring, lopts, NULL)) != -1) {
        switch (opt) {
        case 'f':
            retval->output_format_string = optarg;
            break;
        case 'h':
            retval->unaccounted_header_name = optarg;
            break;
        case 'v':
            retval->verbose = true;
            break;
        case 'u':
            retval->unaccounted_name = optarg;
            break;
        case 'd':
            retval->tag = optarg;
            default_format_string = "%*_%!.%.";
            break;
        case 'M': {
            char *end = optarg;
            retval->max_split = strtol(optarg, &end, 10);
            if (*optarg == '\0' || *end != '\0' || retval->max_split == 0) {
                print_error("split", "Invalid -M argument: \"%s\"", optarg);
                free(retval);
                return NULL;
            }
            if (retval->max_split < 0) // No limit requested
                retval->max_split = LONG_MAX;
            break;
        }
        case 'p': {
            char *end = optarg;
            unsigned long val = strtoul(optarg, &end, 10);
            if (*optarg == '\0' || *end != '\0' || val > 20) {
                print_error("split", "Invalid -p argument: \"%s\"", optarg);
                free(retval);
                return NULL;
            }
            retval->zero_pad = (int) val;
        }
        case 1:
            retval->no_pg = 1;
            break;
        default:
            if (parse_sam_global_opt(opt, optarg, lopts, &retval->ga) == 0) break;
            /* else fall-through */
        case '?':
            usage(stdout);
            free(retval);
            return NULL;
        }
    }

    if (retval->output_format_string == NULL) retval->output_format_string = default_format_string;

    argc -= optind;
    argv += optind;

    if (argc != 1) {
        print_error("split", "Invalid number of arguments: %d", argc);
        usage(stderr);
        free(retval);
        return NULL;
    }

    retval->merged_input_name = argv[0];

    return retval;
}

// Expands a output filename format string
static char* expand_format_string(const char* format_string,
                                  const char* basename, const char* tag_val,
                                  const int file_idx, const int zero_pad,
                                  const htsFormat *format)
{
    kstring_t str = { 0, 0, NULL };
    const char* pointer = format_string;
    const char* next;
    while ((next = strchr(pointer, '%')) != NULL) {
        if (kputsn(pointer, next-pointer, &str) < 0) goto memfail;
        ++next;
        switch (*next) {
            case '%':
                if (kputc('%', &str) < 0) goto memfail;
                break;
            case '*':
                if (kputs(basename, &str) < 0) goto memfail;
                break;
            case '#':
                if (zero_pad == 0) {
                    if (kputl(file_idx, &str) < 0) goto memfail;
                } else {
                    if (ksprintf(&str, "%0*d", zero_pad, file_idx) < 0)
                        goto memfail;
                }
                break;
            case '!':
                if (kputs(tag_val, &str) < 0) goto memfail;
                break;
            case '.':
                // Only really need to cope with sam, bam, cram
                if (format->format != unknown_format) {
                    if (kputs(hts_format_file_extension(format), &str) < 0)
                        goto memfail;
                } else {
                    if (kputs("bam", &str) < 0) goto memfail;
                }
                break;
            case '\0':
                print_error("split", "Trailing %% in filename format string");
                goto fail;
            default:
                // Error is: fprintf(stderr, "bad format string, unknown format specifier\n");
                print_error("split", "Unknown specifier %%%c in filename format string", *next);
                goto fail;
        }
        pointer = next + 1;
    }
    if (kputs(pointer, &str) < 0) goto memfail;
    return ks_release(&str);

 memfail:
    print_error_errno("split", "Couldn't build output filename");
 fail:
    free(str.s);
    return NULL;
}

// Parse the header, count the number of RG tags and return a list of their names
static bool count_RG(sam_hdr_t* hdr, size_t* count, char*** output_name)
{
    char **names = NULL;
    kstring_t id_val = KS_INITIALIZE;
    int i, n_rg = sam_hdr_count_lines(hdr, "RG");

    if (n_rg < 0) {
        print_error("split", "Failed to get @RG IDs");
        *count = 0;
        *output_name = NULL;
        return false;
    }

    if (n_rg == 0) {
        *count = 0;
        *output_name = NULL;
        return true;
    }

    names = calloc(n_rg, sizeof(names[0]));
    if (!names) goto memfail;

    for (i = 0; i < n_rg; i++) {
        if (sam_hdr_find_tag_pos(hdr, "RG", i, "ID", &id_val) < 0) goto memfail;
        names[i] = ks_release(&id_val);
    }

    *count = n_rg;
    *output_name = names;
    return true;

 memfail:
    print_error_errno("split", "Failed to get @RG IDs");
    *count = 0;
    *output_name = NULL;
    ks_free(&id_val);
    free(names);
    return false;
}

static int header_compatible(sam_hdr_t *hdr1, sam_hdr_t *hdr2)
{
    size_t n;
    if (sam_hdr_nref(hdr1) != sam_hdr_nref(hdr2)) {
        print_error("split",
                    "Unaccounted header contains wrong number of references");
        return -1;
    }
    for (n = 0; n < sam_hdr_nref(hdr1); n++) {
        hts_pos_t h1_len = sam_hdr_tid2len(hdr1, n);
        hts_pos_t h2_len = sam_hdr_tid2len(hdr2, n);
        if (h1_len != h2_len) {
            print_error("split",
                        "Unaccounted header reference %zu \"%s\" is not the same length as in the input file",
                        n + 1, sam_hdr_tid2name(hdr2, n));
            return -1;
        }
    }
    return 0;
}

static int grow_output_lists(state_t *state, size_t count) {
    char **new_list = realloc(state->tag_vals, count * sizeof(char *));
    if (!new_list)
        return -1;
    state->tag_vals = new_list;
    new_list = realloc(state->index_file_name, count * sizeof(char *));
    if (!new_list)
        return -1;
    state->index_file_name = new_list;
    new_list = realloc(state->output_file_name, count * sizeof(char *));
    if (!new_list)
        return -1;
    state->output_file_name = new_list;
    samFile **new_file = realloc(state->output_file,
                                 count * sizeof(samFile *));
    if (!new_file)
        return -1;
    state->output_file = new_file;
    sam_hdr_t **new_hdr = realloc(state->output_header,
                                  count * sizeof(sam_hdr_t *));
    if (!new_hdr)
        return -1;
    state->output_header = new_hdr;
    return 0;
}

static khiter_t prep_sam_file(parsed_opts_t *opts, state_t *state,
                              const char *tag, const char *arg_list,
                              int is_rg) {
    char *input_base_name = NULL, *new_file_name = NULL, *tag_val = NULL;
    char *new_idx_fn = NULL;
    sam_hdr_t *new_hdr = NULL;
    samFile *new_sam_file = NULL;

    khiter_t i = kh_get_c2i(state->tag_val_hash, tag);
    if (i != kh_end(state->tag_val_hash)) {
        return i;
    }
    // create new file
    if (grow_output_lists(state, state->output_count + 1) != 0) {
        print_error_errno("split", "Couldn't grow output lists");
        return kh_end(state->tag_val_hash);
    }
    tag_val = strdup(tag);
    if (!tag_val) {
        print_error_errno("split", "Couldn't copy tag value");
        return kh_end(state->tag_val_hash);
    }
    char *dirsep = strrchr(opts->merged_input_name, '/');
    input_base_name = strdup(dirsep? dirsep+1 : opts->merged_input_name);
    if (!input_base_name) {
        print_error_errno("split", "Filename parsing failed");
        goto fail;
    }

    char* extension = strrchr(input_base_name, '.');
    if (extension) *extension = '\0';

    new_file_name = expand_format_string(opts->output_format_string,
                                         input_base_name, tag,
                                         kh_size(state->tag_val_hash),
                                         opts->zero_pad, &opts->ga.out);
    if (!new_file_name) {
        print_error_errno("split", "Filename creation failed");
        goto fail;
    }

    new_hdr = sam_hdr_dup(state->merged_input_header);
    if (!new_hdr) {
        print_error_errno("split", "Duplicating header for file \"%s\" failed", new_file_name);
        goto fail;
    }
    if (!opts->no_pg && sam_hdr_add_pg(new_hdr, "samtools",
                                       "VN", samtools_version(),
                                       arg_list ? "CL": NULL,
                                       arg_list ? arg_list : NULL,
                                       NULL)) {
        print_error_errno("split", "Adding PG line to file \"%s\" failed", new_file_name);
        goto fail;
    }

    if (is_rg) {
        // If here, we've found an RG:Z: tag without a corresponding @RG
        // line in the header.
        if (sam_hdr_remove_lines(new_hdr, "RG", "ID", NULL) != 0) {
            print_error_errno("split",
                              "Failed to remove @RG lines from file \"%s\"",
                              new_file_name);
            goto fail;
        }
        if (sam_hdr_add_line(new_hdr, "RG", "ID", tag_val, NULL) != 0) {
            print_error_errno("split",
                              "Failed to add @RG line to file \"%s\"",
                              new_file_name);
            goto fail;
        }
    }

    char outmode[4] = "w";
    sam_open_mode(outmode + 1, new_file_name, NULL);
    new_sam_file = sam_open_format(new_file_name, outmode, &opts->ga.out);
    if (!new_sam_file) {
        print_error_errno("split", "Opening filename for writing \"%s\" failed", new_file_name);
        goto fail;
    }
    if (state->p.pool)
        hts_set_opt(new_sam_file, HTS_OPT_THREAD_POOL, &state->p);

    if (sam_hdr_write(new_sam_file, new_hdr) != 0) {
        print_error_errno("split", "Couldn't write header to \"%s\"",
                          new_file_name);
        goto fail;
    }

    if (state->write_index) {
        new_idx_fn = auto_index(new_sam_file, new_file_name, new_hdr);
        if (!new_idx_fn) {
            print_error_errno("split", "Creating index file for file \"%s\" failed", new_file_name);
            goto fail;
        }
    }

    int ret = -1;
    i = kh_put_c2i(state->tag_val_hash, tag_val, &ret);
    if (ret < 0) {
        print_error_errno("split", "Adding file \"%s\" failed", new_file_name);
        goto fail;
    }

    kh_val(state->tag_val_hash, i) = state->output_count;
    state->tag_vals[state->output_count] = tag_val;
    state->index_file_name[state->output_count] = new_idx_fn;
    state->output_file_name[state->output_count] = new_file_name;
    state->output_file[state->output_count] = new_sam_file;
    state->output_header[state->output_count] = new_hdr;
    state->output_count++;
    free(input_base_name);
    return i;

 fail:
    free(input_base_name);
    free(new_file_name);
    free(tag_val);
    free(new_idx_fn);
    if (new_hdr)
        sam_hdr_destroy(new_hdr);
    if (new_sam_file)
        sam_close(new_sam_file);
    return kh_end(state->tag_val_hash);
}

// Set the initial state
static state_t* init(parsed_opts_t* opts, const char *arg_list)
{
    state_t* retval = calloc(sizeof(state_t), 1);
    if (!retval) {
        print_error_errno("split", "Initialisation failed");
        return NULL;
    }

    if (opts->ga.nthreads > 0) {
        if (!(retval->p.pool = hts_tpool_init(opts->ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            cleanup_state(retval, false);
            return NULL;
        }
    }

    retval->merged_input_file = sam_open_format(opts->merged_input_name, "r", &opts->ga.in);
    if (!retval->merged_input_file) {
        print_error_errno("split", "Could not open \"%s\"", opts->merged_input_name);
        cleanup_state(retval, false);
        return NULL;
    }
    if (retval->p.pool)
        hts_set_opt(retval->merged_input_file, HTS_OPT_THREAD_POOL, &retval->p);
    retval->merged_input_header = sam_hdr_read(retval->merged_input_file);
    if (retval->merged_input_header == NULL) {
        print_error("split", "Could not read header from \"%s\"", opts->merged_input_name);
        cleanup_state(retval, false);
        return NULL;
    }
    retval->write_index = opts->ga.write_index;

    if (opts->unaccounted_name) {
        if (opts->unaccounted_header_name) {
            samFile* hdr_load = sam_open_format(opts->unaccounted_header_name, "r", &opts->ga.in);
            if (!hdr_load) {
                print_error_errno("split", "Could not open unaccounted header file \"%s\"", opts->unaccounted_header_name);
                cleanup_state(retval, false);
                return NULL;
            }
            retval->unaccounted_header = sam_hdr_read(hdr_load);
            if (retval->unaccounted_header == NULL) {
                print_error("split", "Could not read header from \"%s\"", opts->unaccounted_header_name);
                cleanup_state(retval, false);
                sam_close(hdr_load);
                return NULL;
            }
            sam_close(hdr_load);
            if (header_compatible(retval->merged_input_header,
                                  retval->unaccounted_header) != 0) {
                cleanup_state(retval, false);
                return NULL;
            }
        } else {
            retval->unaccounted_header = sam_hdr_dup(retval->merged_input_header);
            if (!opts->no_pg && sam_hdr_add_pg(retval->unaccounted_header, "samtools",
                                               "VN", samtools_version(),
                                               arg_list ? "CL": NULL,
                                               arg_list ? arg_list : NULL,
                                               NULL)) {
                print_error("split", "Could not rewrite header for \"%s\"", opts->unaccounted_name);
                cleanup_state(retval, false);
                return NULL;
            }
        }

        char outmode[4] = "w";
        sam_open_mode(outmode + 1, opts->unaccounted_name, NULL);
        retval->unaccounted_file = sam_open_format(opts->unaccounted_name, outmode, &opts->ga.out);

        if (retval->unaccounted_file == NULL) {
            print_error_errno("split", "Could not open unaccounted output file \"%s\"", opts->unaccounted_name);
            cleanup_state(retval, false);
            return NULL;
        }
        if (retval->p.pool)
            hts_set_opt(retval->unaccounted_file, HTS_OPT_THREAD_POOL, &retval->p);
    }

    int is_rg = !opts->tag || strcmp(opts->tag, "RG") == 0;
    if (is_rg) {
        if (!count_RG(retval->merged_input_header,
                      &retval->output_count, &retval->tag_vals)) {
            cleanup_state(retval, false);
            return NULL;
        }
        if (opts->verbose)
            fprintf(stderr, "@RG's found %zu\n",retval->output_count);
    } else {
        retval->output_count = 0;
    }

    // Prevent calloc(0, size);
    size_t num = retval->output_count ? retval->output_count : 1;
    retval->index_file_name = (char **)calloc(num, sizeof(char *));
    retval->output_file_name = (char **)calloc(num, sizeof(char *));
    retval->output_file = (samFile**)calloc(num, sizeof(samFile*));
    retval->output_header = (sam_hdr_t**)calloc(num, sizeof(sam_hdr_t*));
    retval->tag_val_hash = kh_init_c2i();
    if (!retval->output_file_name || !retval->output_file || !retval->output_header ||
        !retval->tag_val_hash || !retval->index_file_name) {
        print_error_errno("split", "Could not initialise output file array");
        cleanup_state(retval, false);
        return NULL;
    }
    if (!is_rg)
        return retval;  // Done for this case - outputs will be opened later

    // Adjust max_split if too small for the read-groups listed in the header
    if (opts->max_split < retval->output_count)
        opts->max_split = retval->output_count;

    // Open output files for RGs
    char* dirsep = strrchr(opts->merged_input_name, '/');
    char* input_base_name = strdup(dirsep? dirsep+1 : opts->merged_input_name);
    if (!input_base_name) {
        print_error_errno("split", "Filename manipulation failed");
        cleanup_state(retval, false);
        return NULL;
    }
    char* extension = strrchr(input_base_name, '.');
    if (extension) *extension = '\0';

    size_t i;
    for (i = 0; i < retval->output_count; i++) {
        char* output_filename = NULL;
        char outmode[4] = "w";

        output_filename = expand_format_string(opts->output_format_string,
                                               input_base_name,
                                               retval->tag_vals[i], i,
                                               opts->zero_pad, &opts->ga.out);

        if ( output_filename == NULL ) {
            cleanup_state(retval, false);
            free(input_base_name);
            return NULL;
        }

        retval->output_file_name[i] = output_filename;
        sam_open_mode(outmode + 1, output_filename, NULL);
        retval->output_file[i] = sam_open_format(output_filename, outmode, &opts->ga.out);

        if (retval->output_file[i] == NULL) {
            print_error_errno("split", "Could not open \"%s\"", output_filename);
            cleanup_state(retval, false);
            free(input_base_name);
            return NULL;
        }
        if (retval->p.pool)
            hts_set_opt(retval->output_file[i], HTS_OPT_THREAD_POOL, &retval->p);

        // Record index in hash
        int ret;
        khiter_t iter = kh_put_c2i(retval->tag_val_hash, retval->tag_vals[i], &ret);
        if (ret < 0) {
            print_error_errno("split", "Couldn't add @RG ID to look-up table");
            cleanup_state(retval, false);
            free(input_base_name);
            return NULL;
        }
        kh_val(retval->tag_val_hash,iter) = i;

        // Set and edit header
        retval->output_header[i] = sam_hdr_dup(retval->merged_input_header);
        if (sam_hdr_remove_except(retval->output_header[i], "RG", "ID", retval->tag_vals[i]) ||
            (!opts->no_pg &&
             sam_hdr_add_pg(retval->output_header[i], "samtools",
                            "VN", samtools_version(),
                            arg_list ? "CL": NULL,
                            arg_list ? arg_list : NULL,
                            NULL))) {
            print_error("split", "Could not rewrite header for \"%s\"", output_filename);
            cleanup_state(retval, false);
            free(input_base_name);
            return NULL;
        }
    }

    free(input_base_name);

    return retval;
}

static bool split(state_t* state, parsed_opts_t *opts, char *arg_list)
{
    int is_rg = !opts->tag || strcmp(opts->tag, "RG") == 0;
    if (state->unaccounted_file) {
        if (sam_hdr_write(state->unaccounted_file, state->unaccounted_header) != 0) {
            print_error_errno("split", "Could not write output file header");
            return false;
        }
        if (opts->ga.write_index) {
            state->unaccounted_idx_fn = auto_index(state->unaccounted_file,
                                                   opts->unaccounted_name,
                                                   state->unaccounted_header);
            if (!state->unaccounted_idx_fn) {
                print_error_errno("split",
                                  "Creating index file for file \"%s\" failed",
                                  opts->unaccounted_name);
                return false;
            }
        }
    }
    bam1_t* file_read = bam_init1();
    // Read the first record
    int r;
    if ((r=sam_read1(state->merged_input_file, state->merged_input_header, file_read)) < 0) {
        // Nothing more to read?  Ignore this file
        bam_destroy1(file_read);
        file_read = NULL;
        if (r < -1) {
            print_error("split", "Could not read first input record");
            return false;
        }
    }

    if (is_rg) {
        size_t i;
        for (i = 0; i < state->output_count; i++) {
            if (sam_hdr_write(state->output_file[i], state->output_header[i]) != 0) {
                print_error_errno("split", "Could not write file header to \"%s\"", state->output_file_name[i]);
                goto error;
            }
            if (state->write_index) {
                state->index_file_name[i] = auto_index(state->output_file[i],
                        state->output_file_name[i],
                        state->output_header[i]);
                if (!state->index_file_name[i]) {
                    print_error_errno("split", "Could not create index for file \"%s\"", state->output_file_name[i]);
                    goto error;
                }
            }
        }
    }
    while (file_read != NULL) {
        // Get RG tag from read and look it up in hash to find file to output it to
        uint8_t* tag = bam_aux_get(file_read, is_rg ? "RG" : opts->tag);
        char *val = NULL;
        char number[28];
        khiter_t iter;
        if (tag) {
            switch (*tag) {
            case 'Z': case 'H':
                val = bam_aux2Z(tag);
                break;
            case 'c': case 'C': case 's': case 'S': case 'i': case 'I':
                if (opts->zero_pad == 0) {
                    snprintf(number, sizeof(number), "%"PRId64, bam_aux2i(tag));
                } else {
                    int64_t v = bam_aux2i(tag);
                    snprintf(number, sizeof(number), "%0*"PRId64,
                             v < 0 ? opts->zero_pad + 1 : opts->zero_pad, v);
                }
                val = number;
                break;
            default:
                break;
            }
        }
        if ( val != NULL ) {
            iter = kh_get_c2i(state->tag_val_hash, val);
        } else {
            iter = kh_end(state->tag_val_hash);
        }

        // Check for opts->tag here instead of !is_rg so we open new
        // files if the user specified '-d RG' and we find a RG:Z: value
        // that wasn't listed in the header.  If the '-d' option is
        // not used, we don't open a file to preserve existing behaviour.
        if (opts->tag && val && iter == kh_end(state->tag_val_hash)
            && state->output_count < opts->max_split) {
            // Need to open a new output file
            iter = prep_sam_file(opts, state, val, arg_list, is_rg);
            if (iter == kh_end(state->tag_val_hash)) { // Open failed
                print_error("split",
                            "Could not create output file for tag \"%s:%s\"",
                            opts->tag, bam_aux2Z(tag));
                goto error;

            }
        }

        // Write the read out to correct file
        if (iter != kh_end(state->tag_val_hash)) {
            // if found write to the appropriate untangled bam
            int i = kh_val(state->tag_val_hash,iter);
            if (sam_write1(state->output_file[i], state->output_header[i], file_read) < 0) {
                print_error_errno("split", "Could not write to \"%s\"", state->output_file_name[i]);
                goto error;
            }
        } else {
            // otherwise write to the unaccounted bam if there is one or fail
            if (state->unaccounted_file == NULL) {
                if (tag) {
                    fprintf(stderr, "Read \"%s\" with unaccounted for tag \"%s\".\n", bam_get_qname(file_read), bam_aux2Z(tag));
                } else {
                    fprintf(stderr, "Read \"%s\" has no %s tag.\n",
                            bam_get_qname(file_read), is_rg ? "RG" : opts->tag);
                }
                goto error;
            } else {
                if (sam_write1(state->unaccounted_file, state->unaccounted_header, file_read) < 0) {
                    print_error_errno("split", "Could not write to unaccounted output file");
                    goto error;
                }
            }
        }

        // Replace written read with the next one to process
        if ((r=sam_read1(state->merged_input_file, state->merged_input_header, file_read)) < 0) {
            // Nothing more to read?  Ignore this file in future
            bam_destroy1(file_read);
            file_read = NULL;
            if (r < -1) {
                print_error("split", "Could not read input record");
                return false;
            }
        }
    }

    if (state->write_index) {
        size_t i;
        for (i = 0; i < state->output_count; i++) {
            if (sam_idx_save(state->output_file[i]) < 0) {
                print_error_errno("split", "writing index \"%s\" failed",
                                  state->index_file_name[i]);
                return false;
            }
        }
        if (state->unaccounted_file) {
            if (sam_idx_save(state->unaccounted_file) < 0) {
                print_error_errno("split", "writing index \"%s\" failed",
                                  state->unaccounted_idx_fn);
                return false;
            }
        }
    }

    return true;
error:
    bam_destroy1(file_read);
    return false;
}

static int cleanup_state(state_t* status, bool check_close)
{
    int ret = 0;

    if (!status) return 0;
    if (status->unaccounted_header) sam_hdr_destroy(status->unaccounted_header);
    if (status->unaccounted_file) {
        if (sam_close(status->unaccounted_file) < 0 && check_close) {
            print_error("split", "Error on closing unaccounted file");
            ret = -1;
        }
    }
    sam_close(status->merged_input_file);
    size_t i;
    for (i = 0; i < status->output_count; i++) {
        if (status->output_header && status->output_header[i])
            sam_hdr_destroy(status->output_header[i]);
        if (status->output_file && status->output_file[i]) {
            if (sam_close(status->output_file[i]) < 0 && check_close) {
                print_error("split", "Error on closing output file \"%s\"", status->output_file_name[i]);
                ret = -1;
            }
        }
        if (status->tag_vals) free(status->tag_vals[i]);
        if (status->output_file_name) free(status->output_file_name[i]);
        if (status->index_file_name[i]) free(status->index_file_name[i]);
    }
    if (status->merged_input_header)
        sam_hdr_destroy(status->merged_input_header);
    free(status->output_header);
    free(status->output_file);
    free(status->output_file_name);
    free(status->index_file_name);
    free(status->unaccounted_idx_fn);
    kh_destroy_c2i(status->tag_val_hash);

    free(status->tag_vals);
    if (status->p.pool)
        hts_tpool_destroy(status->p.pool);
    free(status);

    return ret;
}

static void cleanup_opts(parsed_opts_t* opts)
{
    if (!opts) return;
    sam_global_args_free(&opts->ga);
    free(opts);
}

int main_split(int argc, char** argv)
{
    int ret = 1;
    char *arg_list = NULL;
    parsed_opts_t* opts = parse_args(argc, argv);
    if (!opts) goto cleanup_opts;
    if (!opts->no_pg && !(arg_list = stringify_argv(argc+1, argv-1)))
        goto cleanup_opts;
    state_t* status = init(opts, arg_list);
    if (!status) goto cleanup_opts;

    if (!split(status, opts, arg_list)) {
        cleanup_state(status, false);
        goto cleanup_opts;
    }

    ret = cleanup_state(status, true);

cleanup_opts:
    cleanup_opts(opts);
    free(arg_list);

    return ret;
}
