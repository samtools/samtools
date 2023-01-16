/*  bamshuf.c -- collate subcommand.

    Copyright (C) 2012 Broad Institute.
    Copyright (C) 2013, 2015-2019,2023 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#endif
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/ksort.h"
#include "samtools.h"
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include "htslib/khash.h"

#define DEF_CLEVEL 1

static inline unsigned hash_Wang(unsigned key)
{
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}

static inline unsigned hash_X31_Wang(const char *s)
{
    unsigned h = *s;
    if (h) {
        for (++s ; *s; ++s) h = (h << 5) - h + *s;
        return hash_Wang(h);
    } else return 0;
}

typedef struct {
    unsigned key;
    bam1_t *b;
} elem_t;

static inline int elem_lt(elem_t x, elem_t y)
{
    if (x.key < y.key) return 1;
    if (x.key == y.key) {
        int t;
        t = strcmp(bam_get_qname(x.b), bam_get_qname(y.b));
        if (t < 0) return 1;
        return (t == 0 && ((x.b->core.flag>>6&3) < (y.b->core.flag>>6&3)));
    } else return 0;
}

KSORT_INIT(bamshuf, elem_t, elem_lt)


typedef struct {
    int written;
    bam1_t *b;
} bam_item_t;

typedef struct {
    bam1_t *bam_pool;
    bam_item_t *items;
    size_t size;
    size_t index;
} bam_list_t;

typedef struct {
    bam_item_t *bi;
} store_item_t;

KHASH_MAP_INIT_STR(bam_store, store_item_t)


static bam_item_t *store_bam(bam_list_t *list) {
    size_t old_index = list->index;

    list->items[list->index++].written = 0;

    if (list->index >= list->size)
        list->index = 0;

    return &list->items[old_index];
}


static int write_bam_needed(bam_list_t *list) {
    return !list->items[list->index].written;
}


static void mark_bam_as_written(bam_list_t *list) {
    list->items[list->index].written = 1;
}


static int create_bam_list(bam_list_t *list, size_t max_size) {
    size_t i;

    list->size = list->index = 0;
    list->items    = NULL;
    list->bam_pool = NULL;

    if ((list->items = malloc(max_size * sizeof(bam_item_t))) == NULL) {
        return 1;
    }

    if ((list->bam_pool = calloc(max_size, sizeof(bam1_t))) == NULL) {
        return 1;
    }

    for (i = 0; i < max_size; i++) {
        list->items[i].b = &list->bam_pool[i];
        list->items[i].written = 1;
    }

    list->size  = max_size;
    list->index = 0;

    return 0;
}


static void destroy_bam_list(bam_list_t *list) {
    size_t i;

    for (i = 0; i < list->size; i++) {
        free(list->bam_pool[i].data);
    }

    free(list->bam_pool);
    free(list->items);
}


static inline int write_to_bin_file(bam1_t *bam, int64_t *count, samFile **bin_files, char **names, sam_hdr_t *header, int files) {
    uint32_t x;

    x = hash_X31_Wang(bam_get_qname(bam)) % files;

    if (sam_write1(bin_files[x], header, bam) < 0) {
        print_error_errno("collate", "Couldn't write to intermediate file \"%s\"", names[x]);
        return 1;
    }

    ++count[x];

    return 0;
}


static int bamshuf(const char *fn, int n_files, const char *pre, int clevel,
                   int is_stdout, const char *output_file, int fast, int store_max, sam_global_args *ga, char *arg_list, int no_pg)
{
    samFile *fp, *fpw = NULL, **fpt = NULL;
    char **fnt = NULL, modew[8];
    bam1_t *b = NULL;
    int i, counter, l, r;
    sam_hdr_t *h = NULL;
    int64_t j, max_cnt = 0, *cnt = NULL;
    elem_t *a = NULL;
    htsThreadPool p = {NULL, 0};

    if (ga->nthreads > 0) {
        if (!(p.pool = hts_tpool_init(ga->nthreads))) {
            print_error_errno("collate", "Error creating thread pool\n");
            return 1;
        }
    }

    // Read input, distribute reads pseudo-randomly into n_files temporary
    // files.
    fp = sam_open_format(fn ? fn : "-", "r", &ga->in);
    if (fp == NULL) {
        print_error_errno("collate", "Cannot open input file \"%s\"", fn);
        return 1;
    }
    if (p.pool) hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);

    h = sam_hdr_read(fp);
    if (h == NULL) {
        fprintf(stderr, "Couldn't read header for '%s'\n", fn);
        goto fail;
    }

    if ((-1 == sam_hdr_update_hd(h, "SO", "unsorted", "GO", "query"))
     && (-1 == sam_hdr_add_line(h, "HD", "VN", SAM_FORMAT_VERSION, "SO", "unsorted", "GO", "query", NULL))
     ) {
        print_error("collate", "failed to update HD line\n");
        goto fail;
    }

    // open final output file
    l = strlen(pre);

    sprintf(modew, "wb%d", (clevel >= 0 && clevel <= 9)? clevel : DEF_CLEVEL);

    if (!is_stdout && !output_file) { // output to a file (name based on prefix)
        char *fnw = (char*)calloc(l + 5, 1);
        if (!fnw) goto mem_fail;
        if (ga->out.format == unknown_format)
            sprintf(fnw, "%s.bam", pre); // "wb" above makes BAM the default
        else
            sprintf(fnw, "%s.%s", pre,  hts_format_file_extension(&ga->out));
        fpw = sam_open_format(fnw, modew, &ga->out);
        free(fnw);
    } else if (output_file) { // output to a given file
        modew[0] = 'w'; modew[1] = '\0';
        sam_open_mode(modew + 1, output_file, NULL);
        j = strlen(modew);
        snprintf(modew + j, sizeof(modew) - j, "%d",
                 (clevel >= 0 && clevel <= 9)? clevel : DEF_CLEVEL);
        fpw = sam_open_format(output_file, modew, &ga->out);
    } else fpw = sam_open_format("-", modew, &ga->out); // output to stdout
    if (fpw == NULL) {
        if (is_stdout) print_error_errno("collate", "Cannot open standard output");
        else print_error_errno("collate", "Cannot open output file \"%s.bam\"", pre);
        goto fail;
    }
    if (p.pool) hts_set_opt(fpw, HTS_OPT_THREAD_POOL, &p);

    if (!no_pg && sam_hdr_add_pg(h, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL)) {
        print_error("collate", "failed to add PG line to header of \"%s\"", output_file);
        goto fail;
    }

    if (sam_hdr_write(fpw, h) < 0) {
        print_error_errno("collate", "Couldn't write header");
        goto fail;
    }

    fnt = (char**)calloc(n_files, sizeof(char*));
    if (!fnt) goto mem_fail;
    fpt = (samFile**)calloc(n_files, sizeof(samFile*));
    if (!fpt) goto mem_fail;
    cnt = (int64_t*)calloc(n_files, 8);
    if (!cnt) goto mem_fail;

    for (i = counter = 0; i < n_files; ++i) {
        fnt[i] = (char*)calloc(l + 20, 1);
        if (!fnt[i]) goto mem_fail;
        do {
            sprintf(fnt[i], "%s.%04d.bam", pre, counter++);
            fpt[i] = sam_open(fnt[i], "wxb1");
        } while (!fpt[i] && errno == EEXIST);
        if (fpt[i] == NULL) {
            print_error_errno("collate", "Cannot open intermediate file \"%s\"", fnt[i]);
            goto fail;
        }
        if (p.pool) hts_set_opt(fpt[i], HTS_OPT_THREAD_POOL, &p);
        if (sam_hdr_write(fpt[i], h) < 0) {
            print_error_errno("collate", "Couldn't write header to intermediate file \"%s\"", fnt[i]);
            goto fail;
        }
    }

    if (fast) {
        khash_t(bam_store) *stored = kh_init(bam_store);
        khiter_t itr;
        bam_list_t list;
        int err = 0;
        if (!stored) goto mem_fail;

        if (store_max < 2) store_max = 2;

        if (create_bam_list(&list, store_max)) {
            fprintf(stderr, "[collate[ ERROR: unable to create bam list.\n");
            err = 1;
            goto fast_fail;
        }

        while ((r = sam_read1(fp, h, list.items[list.index].b)) >= 0) {
            int ret;
            bam1_t *b = list.items[list.index].b;
            int readflag = b->core.flag & (BAM_FREAD1 | BAM_FREAD2);

            // strictly paired reads only
            if (!(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) && (readflag == BAM_FREAD1 || readflag == BAM_FREAD2)) {

                itr = kh_get(bam_store, stored, bam_get_qname(b));

                if (itr == kh_end(stored)) {
                    // new read
                    itr = kh_put(bam_store, stored, bam_get_qname(b), &ret);

                    if (ret > 0) { // okay to go ahead store it
                        kh_value(stored, itr).bi = store_bam(&list);

                        // see if the next one on the list needs to be written out
                        if (write_bam_needed(&list)) {
                            if (write_to_bin_file(list.items[list.index].b, cnt, fpt, fnt, h, n_files) < 0) {
                                fprintf(stderr, "[collate] ERROR: could not write line.\n");
                                err = 1;
                                goto fast_fail;
                            } else {
                                mark_bam_as_written(&list);

                                itr = kh_get(bam_store, stored, bam_get_qname(list.items[list.index].b));

                                if (itr != kh_end(stored)) {
                                    kh_del(bam_store, stored, itr);
                                } else {
                                    fprintf(stderr, "[collate] ERROR: stored value not in hash.\n");
                                    err = 1;
                                    goto fast_fail;
                                }
                            }
                        }
                    } else if (ret == 0) {
                        fprintf(stderr, "[collate] ERROR: value already in hash.\n");
                        err = 1;
                        goto fast_fail;
                    } else {
                        fprintf(stderr, "[collate] ERROR: unable to store in hash.\n");
                        err = 1;
                        goto fast_fail;
                    }
                } else { // we have a match
                    // write out the reads in R1 R2 order
                    bam1_t *r1, *r2;

                    if (b->core.flag & BAM_FREAD1) {
                        r1 = b;
                        r2 = kh_value(stored, itr).bi->b;
                    } else {
                        r1 = kh_value(stored, itr).bi->b;
                        r2 = b;
                    }

                    if (sam_write1(fpw, h, r1) < 0) {
                        fprintf(stderr, "[collate] ERROR: could not write r1 alignment.\n");
                        err = 1;
                        goto fast_fail;
                    }

                    if (sam_write1(fpw, h, r2) < 0) {
                        fprintf(stderr, "[collate] ERROR: could not write r2 alignment.\n");
                        err = 1;
                        goto fast_fail;
                    }

                    mark_bam_as_written(&list);

                    // remove stored read
                    kh_value(stored, itr).bi->written = 1;
                    kh_del(bam_store, stored, itr);
                }
            }
        }

        for (list.index = 0; list.index < list.size; list.index++) {
            if (write_bam_needed(&list)) {
                bam1_t *b = list.items[list.index].b;

                if (write_to_bin_file(b, cnt, fpt, fnt, h, n_files)) {
                    err = 1;
                    goto fast_fail;
                } else {
                    itr = kh_get(bam_store, stored, bam_get_qname(b));
                    kh_del(bam_store, stored, itr);
                }
            }
        }

 fast_fail:
        if (err) {
            for (itr = kh_begin(stored); itr != kh_end(stored); ++itr) {
                if (kh_exist(stored, itr)) {
                    kh_del(bam_store, stored, itr);
                }
            }

            kh_destroy(bam_store, stored);
            destroy_bam_list(&list);
            goto fail;
        } else {
            kh_destroy(bam_store, stored);
            destroy_bam_list(&list);
        }

    } else {
        b = bam_init1();
        if (!b) goto mem_fail;

        while ((r = sam_read1(fp, h, b)) >= 0) {
            if (write_to_bin_file(b, cnt, fpt, fnt, h, n_files)) {
                bam_destroy1(b);
                goto fail;
            }
        }

        bam_destroy1(b);
    }

    if (r < -1) {
        fprintf(stderr, "Error reading input file\n");
        goto fail;
    }
    for (i = 0; i < n_files; ++i) {
        // Close split output
        r = sam_close(fpt[i]);
        fpt[i] = NULL;
        if (r < 0) {
            fprintf(stderr, "Error on closing '%s'\n", fnt[i]);
            return 1;
        }

        // Find biggest count
        if (max_cnt < cnt[i]) max_cnt = cnt[i];
    }
    free(fpt);
    fpt = NULL;
    sam_close(fp);
    fp = NULL;

    // merge
    a = malloc(max_cnt * sizeof(elem_t));
    if (!a) goto mem_fail;
    for (j = 0; j < max_cnt; ++j) {
        a[j].b = bam_init1();
        if (!a[j].b) { max_cnt = j; goto mem_fail; }
    }

    for (i = 0; i < n_files; ++i) {
        int64_t c = cnt[i];
        fp = sam_open_format(fnt[i], "r", &ga->in);
        if (NULL == fp) {
            print_error_errno("collate", "Couldn't open \"%s\"", fnt[i]);
            goto fail;
        }
        if (p.pool) hts_set_opt(fp, HTS_OPT_THREAD_POOL, &p);
        sam_hdr_destroy(sam_hdr_read(fp)); // Skip over header

        // Slurp in one of the split files
        for (j = 0; j < c; ++j) {
            if (sam_read1(fp, h, a[j].b) < 0) {
                fprintf(stderr, "Error reading '%s'\n", fnt[i]);
                goto fail;
            }
            a[j].key = hash_X31_Wang(bam_get_qname(a[j].b));
        }
        sam_close(fp);
        unlink(fnt[i]);
        free(fnt[i]);
        fnt[i] = NULL;

        ks_introsort(bamshuf, c, a); // Shuffle all the reads

        // Write them out again
        for (j = 0; j < c; ++j) {
            if (sam_write1(fpw, h, a[j].b) < 0) {
                print_error_errno("collate", "Error writing to output");
                goto fail;
            }
        }
    }

    sam_hdr_destroy(h);
    for (j = 0; j < max_cnt; ++j) bam_destroy1(a[j].b);
    free(a); free(fnt); free(cnt);
    sam_global_args_free(ga);
    if (sam_close(fpw) < 0) {
        fprintf(stderr, "Error on closing output\n");
        return 1;
    }

    if (p.pool) hts_tpool_destroy(p.pool);
    return 0;

 mem_fail:
    fprintf(stderr, "Out of memory\n");

 fail:
    if (fp) sam_close(fp);
    if (fpw) sam_close(fpw);
    if (h) sam_hdr_destroy(h);
    for (i = 0; i < n_files; ++i) {
        if (fnt) free(fnt[i]);
        if (fpt && fpt[i]) sam_close(fpt[i]);
    }
    if (a) {
        for (j = 0; j < max_cnt; ++j) bam_destroy1(a[j].b);
        free(a);
    }
    free(fnt);
    free(fpt);
    free(cnt);
    if (p.pool) hts_tpool_destroy(p.pool);
    sam_global_args_free(ga);
    return 1;
}

static int usage(FILE *fp, int n_files, int reads_store) {
    fprintf(fp,
            "Usage: samtools collate [options...] <in.bam> [<prefix>]\n\n"
            "Options:\n"
            "      -O       Output to stdout\n"
            "      -o       Output file name (use prefix if not set)\n"
            "      -u       Uncompressed BAM output\n"
            "      -f       Fast (only primary alignments)\n"
            "      -r       Working reads stored (with -f) [%d]\n" // reads_store
            "      -l INT   Compression level [%d]\n" // DEF_CLEVEL
            "      -n INT   Number of temporary files [%d]\n" // n_files
            "      -T PREFIX\n"
            "               Write tempory files to PREFIX.nnnn.bam\n"
            "      --no-PG  do not add a PG line\n",
            reads_store, DEF_CLEVEL, n_files);

    sam_global_opt_help(fp, "-....@-.");
    fprintf(fp,
            "  <prefix> is required unless the -o or -O options are used.\n");

    return 1;
}

char *generate_prefix(const char *out_fn) {
    char *prefix;
    unsigned int pid = getpid();

    if (out_fn && !(*out_fn == '-' && out_fn[1] == '\0')) {
        // <out_fn>.<collate><pid>.<nnnn>.<bam>
        size_t plen = strlen(out_fn) + 50;
        if (!(prefix = malloc(plen))) {
            perror("collate");
            return NULL;
        }
        snprintf(prefix, plen, "%s.collate%x", out_fn, pid);
        return prefix;
    }

#ifdef _WIN32
#  define PREFIX_LEN (MAX_PATH + 16)
    DWORD ret;
    prefix = calloc(PREFIX_LEN, sizeof(*prefix));
    if (!prefix) {
        perror("collate");
        return NULL;
    }
    ret = GetTempPathA(MAX_PATH, prefix);
    if (ret > MAX_PATH || ret == 0) {
        fprintf(stderr,
                "[E::collate] Couldn't get path for temporary files.\n");
        free(prefix);
        return NULL;
    }
    snprintf(prefix + ret, PREFIX_LEN - ret, "\\%x", pid);
    return prefix;
#else
    char *tmp_env = getenv("TMPDIR");
    if (!tmp_env)
        tmp_env = "/tmp";

    size_t prefix_len = strlen(tmp_env)+20;
    prefix = malloc(prefix_len);
    if (!prefix) {
        perror("collate");
        return NULL;
    }
    snprintf(prefix, prefix_len, "%s/collate%x", tmp_env, pid);

    return prefix;
#endif
}

int main_bamshuf(int argc, char *argv[])
{
    int c, n_files = 64, clevel = DEF_CLEVEL, is_stdout = 0, is_un = 0, fast_coll = 0, reads_store = 10000, ret, pre_mem = 0, no_pg = 0;
    const char *output_file = NULL;
    char *prefix = NULL, *arg_list = NULL;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0, '@'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "n:l:uOo:@:fr:T:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'n': n_files = atoi(optarg); break;
        case 'l': clevel = atoi(optarg); break;
        case 'u': is_un = 1; break;
        case 'O': is_stdout = 1; break;
        case 'o': output_file = optarg; break;
        case 'f': fast_coll = 1; break;
        case 'r': reads_store = atoi(optarg); break;
        case 'T': prefix = optarg; break;
        case 1: no_pg = 1; break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                  /* else fall-through */
        case '?': return usage(stderr, n_files, reads_store);
        }
    }
    if (is_un) clevel = 0;
    if (argc >= optind + 2) prefix = argv[optind+1];
    if (argc == optind) {
        if (argc > 1 || !isatty(STDIN_FILENO))
            fprintf(stderr, "collate: no input filename specified.\n");
        return usage(argc > 1 || !isatty(STDIN_FILENO) ? stderr : stdout,
                     n_files, reads_store);
    }
    if (!(prefix || is_stdout || output_file))
        return usage(stderr, n_files, reads_store);
    if (is_stdout && output_file) {
        fprintf(stderr, "collate: -o and -O options cannot be used together.\n");
        return usage(stderr, n_files, reads_store);
    }
    if (!prefix) {
        prefix = generate_prefix(output_file);
        pre_mem = 1;
    }

    if (!prefix) return EXIT_FAILURE;

    if (!no_pg && !(arg_list = stringify_argv(argc+1, argv-1))) {
        print_error("collate", "failed to create arg_list");
        return 1;
    }

    ret = bamshuf(argv[optind], n_files, prefix, clevel, is_stdout,
                   output_file, fast_coll, reads_store, &ga, arg_list, no_pg);

    if (pre_mem) free(prefix);
    free(arg_list);

    return ret;
}
