/*  bamshuf.c -- collate subcommand.

    Copyright (C) 2012 Broad Institute.
    Copyright (C) 2013 Genome Research Ltd.

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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/ksort.h"
#include "samtools.h"
#include "sam_opts.h"

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

static int bamshuf(const char *fn, int n_files, const char *pre, int clevel,
                   int is_stdout, sam_global_args *ga)
{
    samFile *fp, *fpw, **fpt;
    char **fnt, modew[8];
    bam1_t *b;
    int i, l;
    bam_hdr_t *h;
    int64_t *cnt;

    // split
    fp = sam_open_format(fn, "r", &ga->in);
    if (fp == NULL) {
        print_error_errno("collate", "Cannot open input file \"%s\"", fn);
        return 1;
    }

    h = sam_hdr_read(fp);
    if (h == NULL) {
        fprintf(stderr, "Couldn't read header for '%s'\n", fn);
        return 1;
    }
    fnt = (char**)calloc(n_files, sizeof(char*));
    fpt = (samFile**)calloc(n_files, sizeof(samFile*));
    cnt = (int64_t*)calloc(n_files, 8);
    l = strlen(pre);

    for (i = 0; i < n_files; ++i) {
        fnt[i] = (char*)calloc(l + 10, 1);
        sprintf(fnt[i], "%s.%.4d.bam", pre, i);
        fpt[i] = sam_open(fnt[i], "wb1");
        if (fpt[i] == NULL) {
            print_error_errno("collate", "Cannot open intermediate file \"%s\"", fnt[i]);
            return 1;
        }
        sam_hdr_write(fpt[i], h);
    }
    b = bam_init1();
    while (sam_read1(fp, h, b) >= 0) {
        uint32_t x;
        x = hash_X31_Wang(bam_get_qname(b)) % n_files;
        sam_write1(fpt[x], h, b);
        ++cnt[x];
    }
    bam_destroy1(b);
    for (i = 0; i < n_files; ++i) sam_close(fpt[i]);
    free(fpt);
    sam_close(fp);

    // merge
    sprintf(modew, "wb%d", (clevel >= 0 && clevel <= 9)? clevel : DEF_CLEVEL);
    if (!is_stdout) { // output to a file
        char *fnw = (char*)calloc(l + 5, 1);
        if (ga->out.format == unknown_format)
            sprintf(fnw, "%s.bam", pre); // "wb" above makes BAM the default
        else
            sprintf(fnw, "%s.%s", pre,  hts_format_file_extension(&ga->out));
        fpw = sam_open_format(fnw, modew, &ga->out);
        free(fnw);
    } else fpw = sam_open_format("-", modew, &ga->out); // output to stdout
    if (fpw == NULL) {
        if (is_stdout) print_error_errno("collate", "Cannot open standard output");
        else print_error_errno("collate", "Cannot open output file \"%s.bam\"", pre);
        return 1;
    }

    sam_hdr_write(fpw, h);
    for (i = 0; i < n_files; ++i) {
        int64_t j, c = cnt[i];
        elem_t *a;
        fp = sam_open_format(fnt[i], "r", &ga->in);
        bam_hdr_destroy(sam_hdr_read(fp));
        a = (elem_t*)calloc(c, sizeof(elem_t));
        for (j = 0; j < c; ++j) {
            a[j].b = bam_init1();
            sam_read1(fp, h, a[j].b);
            a[j].key = hash_X31_Wang(bam_get_qname(a[j].b));
        }
        sam_close(fp);
        unlink(fnt[i]);
        free(fnt[i]);
        ks_introsort(bamshuf, c, a);
        for (j = 0; j < c; ++j) {
            sam_write1(fpw, h, a[j].b);
            bam_destroy1(a[j].b);
        }
        free(a);
    }
    sam_close(fpw);
    bam_hdr_destroy(h);
    free(fnt); free(cnt);
    sam_global_args_free(ga);

    return 0;
}

static int usage(FILE *fp, int n_files) {
    fprintf(fp,
            "Usage:   samtools collate [-Ou] [-n nFiles] [-c cLevel] <in.bam> <out.prefix>\n\n"
            "Options:\n"
            "      -O       output to stdout\n"
            "      -u       uncompressed BAM output\n"
            "      -l INT   compression level [%d]\n" // DEF_CLEVEL
            "      -n INT   number of temporary files [%d]\n", // n_files
            DEF_CLEVEL, n_files);

    sam_global_opt_help(fp, "-....");

    return 1;
}

int main_bamshuf(int argc, char *argv[])
{
    int c, n_files = 64, clevel = DEF_CLEVEL, is_stdout = 0, is_un = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "n:l:uO", lopts, NULL)) >= 0) {
        switch (c) {
        case 'n': n_files = atoi(optarg); break;
        case 'l': clevel = atoi(optarg); break;
        case 'u': is_un = 1; break;
        case 'O': is_stdout = 1; break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                  /* else fall-through */
        case '?': return usage(stderr, n_files);
        }
    }
    if (is_un) clevel = 0;
    if (optind + 2 > argc)
        return usage(stderr, n_files);

    return bamshuf(argv[optind], n_files, argv[optind+1], clevel, is_stdout, &ga);
}
