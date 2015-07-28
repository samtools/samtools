/*  bamshuf.c -- bamshuf subcommand.

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
#include "htslib/bgzf.h"
#include "htslib/ksort.h"
#include "samtools.h"

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

static int bamshuf(const char *fn, int n_files, const char *pre, int clevel, int is_stdout)
{
    BGZF *fp, *fpw = NULL, **fpt = NULL;
    char **fnt = NULL, modew[8];
    bam1_t *b = NULL;
    int i, l, r;
    bam_hdr_t *h = NULL;
    int64_t *cnt = NULL;
    elem_t *a = NULL;
    int64_t a_size = 0, a_first = 0, a_last = 0;

    // Read input, distribute reads pseudo-randomly into n_files temporary
    // files.
    fp = strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
    if (fp == NULL) {
        print_error_errno("Cannot open input file \"%s\"", fn);
        return 1;
    }

    h = bam_hdr_read(fp);
    if (h == NULL) {
        fprintf(stderr, "Couldn't read header for '%s'\n", fn);
        goto fail;
    }
    fnt = (char**)calloc(n_files, sizeof(char*));
    if (!fnt) goto mem_fail;
    fpt = (BGZF**)calloc(n_files, sizeof(BGZF*));
    if (!fpt) goto mem_fail;
    cnt = (int64_t*)calloc(n_files, 8);
    if (!cnt) goto mem_fail;

    l = strlen(pre);
    for (i = 0; i < n_files; ++i) {
        fnt[i] = (char*)calloc(l + 10, 1);
        if (!fnt[i]) goto mem_fail;
        sprintf(fnt[i], "%s.%.4d.bam", pre, i);
        fpt[i] = bgzf_open(fnt[i], "w1");
        if (fpt[i] == NULL) {
            print_error_errno("Cannot open intermediate file \"%s\"", fnt[i]);
            goto fail;
        }
        if (bam_hdr_write(fpt[i], h) < 0) {
            fprintf(stderr, "Couldn't write header for '%s'\n", fnt[i]);
            goto fail;
        }
    }
    b = bam_init1();
    if (!b) goto mem_fail;
    while ((r = bam_read1(fp, b)) >= 0) {
        uint32_t x;
        x = hash_X31_Wang(bam_get_qname(b)) % n_files;
        if (bam_write1(fpt[x], b) < 0) {
            fprintf(stderr, "Couldn't write to '%s'\n", fnt[x]);
            goto fail;
        }
        ++cnt[x];
    }
    if (r < -1) {
        fprintf(stderr, "Error reading input file\n");
        goto fail;
    } 
    bam_destroy1(b); b = NULL;
    for (i = 0; i < n_files; ++i) {
        // Close split output
        r = bgzf_close(fpt[i]);
        fpt[i] = NULL;
        if (r < 0) {
            fprintf(stderr, "Error on closing '%s'\n", fnt[i]);
            return 1;
        }

        // Find biggest count
        if (a_size < cnt[i]) a_size = cnt[i];
    }
    free(fpt);
    fpt = NULL;
    bgzf_close(fp);
    fp = NULL;
    // merge
    sprintf(modew, "w%d", (clevel >= 0 && clevel <= 9)? clevel : DEF_CLEVEL);
    if (!is_stdout) { // output to a file
        char *fnw = (char*)calloc(l + 5, 1);
        if (!fnw) goto mem_fail;
        sprintf(fnw, "%s.bam", pre);
        fpw = bgzf_open(fnw, modew);
        free(fnw);
    } else fpw = bgzf_dopen(fileno(stdout), modew); // output to stdout
    if (fpw == NULL) {
        if (is_stdout) print_error_errno("Cannot open standard output");
        else print_error_errno("Cannot open output file \"%s.bam\"", pre);
        goto fail;
    }

    if (bam_hdr_write(fpw, h) < 0) {
        fprintf(stderr, "Couldn't write bam header\n");
        goto fail;
    }
    bam_hdr_destroy(h);
    h = NULL;

    a = malloc(a_size * sizeof(elem_t));
    if (!a) goto mem_fail;

    for (i = 0; i < n_files; ++i) {
        int64_t j, c = cnt[i];
        fp = bgzf_open(fnt[i], "r");
        if (NULL == fp) {
            fprintf(stderr, "Couldn't open '%s'\n", fnt[i]);
            goto fail;
        }
        bam_hdr_destroy(bam_hdr_read(fp)); // Skip over header

        a_first = a_last = 0; // Track which bit of 'a' contains bam structs

        // Slurp in one of the split files
        for (j = 0; j < c; ++j) {
            a[j].b = bam_init1();
            if (!a[j].b) goto mem_fail;
            a_last++;
            if (bam_read1(fp, a[j].b) < 0) {
                fprintf(stderr, "Error reading '%s'\n", fnt[i]);
                goto fail;
            }
            a[j].key = hash_X31_Wang(bam_get_qname(a[j].b));
        }
        bgzf_close(fp);
        unlink(fnt[i]);
        free(fnt[i]);
        fnt[i] = NULL;

        ks_introsort(bamshuf, c, a); // Shuffle all the reads

        // Write them out again
        for (j = 0; j < c; ++j) {
            if (bam_write1(fpw, a[j].b) < 0) {
                fprintf(stderr, "Error writing to output\n");
                goto fail;
            }
            bam_destroy1(a[j].b);
            a_first++;
        }
    }
    
    free(a); free(fnt); free(cnt);
    if (bgzf_close(fpw) < 0) {
        fprintf(stderr, "Error on closing output\n");
        return 1;
    }
    return 0;

 mem_fail:
    fprintf(stderr, "Out of memory\n");

 fail:
    if (fp) bgzf_close(fp);
    if (fpw) bgzf_close(fpw);
    if (h) bam_hdr_destroy(h);
    if (b) bam_destroy1(b);
    for (i = 0; i < n_files; ++i) {
        if (fnt) free(fnt[i]);
        if (fpt && fpt[i]) bgzf_close(fpt[i]);
    }
    if (a) {
        int64_t j;
        for (j = a_first; j < a_last; j++) {
            if (a[j].b) bam_destroy1(a[j].b);
        }
        free(a);
    }
    free(fnt);
    free(fpt);
    free(cnt);
    return 1;
}

int main_bamshuf(int argc, char *argv[])
{
    int c, n_files = 64, clevel = DEF_CLEVEL, is_stdout = 0, is_un = 0;
    while ((c = getopt(argc, argv, "n:l:uO")) >= 0) {
        switch (c) {
        case 'n': n_files = atoi(optarg); break;
        case 'l': clevel = atoi(optarg); break;
        case 'u': is_un = 1; break;
        case 'O': is_stdout = 1; break;
        }
    }
    if (is_un) clevel = 0;
    if (optind + 2 > argc) {
        fprintf(stderr,
"Usage:   samtools bamshuf [-Ou] [-n nFiles] [-c cLevel] <in.bam> <out.prefix>\n\n"
"Options: -O      output to stdout\n"
"         -u      uncompressed BAM output\n"
"         -l INT  compression level [%d]\n" // DEF_CLEVEL
"         -n INT  number of temporary files [%d]\n", // n_files
                DEF_CLEVEL, n_files);
        return 1;
    }
    return bamshuf(argv[optind], n_files, argv[optind+1], clevel, is_stdout);
}
