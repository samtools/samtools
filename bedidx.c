/*  bedidx.c -- BED file indexing.

    Copyright (C) 2011 Broad Institute.
    Copyright (C) 2014 Genome Research Ltd.

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

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <zlib.h>

#ifdef _WIN32
#define drand48() ((double)rand() / RAND_MAX)
#endif

#include "htslib/ksort.h"
KSORT_INIT_GENERIC(uint64_t)

#include "htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 8192)

typedef struct {
    int n, m;
    uint64_t *a;
    int *idx;
} bed_reglist_t;

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(reg, bed_reglist_t)

#define LIDX_SHIFT 13

typedef kh_reg_t reghash_t;

void bed_destroy(void *_h);


int *bed_index_core(int n, uint64_t *a, int *n_idx)
{
    int i, j, m, *idx;
    m = *n_idx = 0; idx = 0;
    for (i = 0; i < n; ++i) {
        int beg, end;
        beg = a[i]>>32 >> LIDX_SHIFT; end = ((uint32_t)a[i]) >> LIDX_SHIFT;
        if (m < end + 1) {
            int oldm = m;
            m = end + 1;
            kroundup32(m);
            idx = realloc(idx, m * sizeof(int));
            for (j = oldm; j < m; ++j) idx[j] = -1;
        }
        if (beg == end) {
            if (idx[beg] < 0) idx[beg] = i;
        } else {
            for (j = beg; j <= end; ++j)
                if (idx[j] < 0) idx[j] = i;
        }
        *n_idx = end + 1;
    }
    return idx;
}

void bed_index(void *_h)
{
    reghash_t *h = (reghash_t*)_h;
    khint_t k;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            bed_reglist_t *p = &kh_val(h, k);
            if (p->idx) free(p->idx);
            ks_introsort(uint64_t, p->n, p->a);
            p->idx = bed_index_core(p->n, p->a, &p->m);
        }
    }
}

int bed_overlap_core(const bed_reglist_t *p, int beg, int end)
{
    int i, min_off;
    if (p->n == 0) return 0;
    min_off = (beg>>LIDX_SHIFT >= p->n)? p->idx[p->n-1] : p->idx[beg>>LIDX_SHIFT];
    if (min_off < 0) { // TODO: this block can be improved, but speed should not matter too much here
        int n = beg>>LIDX_SHIFT;
        if (n > p->n) n = p->n;
        for (i = n - 1; i >= 0; --i)
            if (p->idx[i] >= 0) break;
        min_off = i >= 0? p->idx[i] : 0;
    }
    for (i = min_off; i < p->n; ++i) {
        if ((int)(p->a[i]>>32) >= end) break; // out of range; no need to proceed
        if ((int32_t)p->a[i] > beg && (int32_t)(p->a[i]>>32) < end)
            return 1; // find the overlap; return
    }
    return 0;
}

int bed_overlap(const void *_h, const char *chr, int beg, int end)
{
    const reghash_t *h = (const reghash_t*)_h;
    khint_t k;
    if (!h) return 0;
    k = kh_get(reg, h, chr);
    if (k == kh_end(h)) return 0;
    return bed_overlap_core(&kh_val(h, k), beg, end);
}

/* "BED" file reader, which actually reads two different formats.

   BED files contain between three and nine fields per line, of which
   only the first three (reference, start, end) are of interest to us.
   BED counts positions from base 0, and the end is the base after the
   region of interest.  While not properly documented in the specification,
   it is also possible to have 'browser' and 'track' lines in BED files that
   do not follow the standard format and should be ignored.  Examination
   of the BED file reading code in
   http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git shows that BED
   files can also have comment lines starting with '#', leading whitespace
   is stripped, and that fields are separated by one or more consecutive
   whitespace characters.

   The alternative format was originally for reading positions in VCF
   format.  This expects two columns, which indicate the reference and
   a position.  The position corresponds to a single base, and unlike
   BED counts from 1.

   Which format is in use is determined based on whether one or two
   numbers can be decoded on the line.  As this choice is made line-by-line
   in this implementation, it is possible (but probably a bad idea) to mix
   both formats in the same file.  If trying to read a VCF file by this
   method, it would be important to ensure that the third column (ID) does
   not contain any entries that start with a digit, to avoid the line
   erroneously being parsed as a BED file entry.

   The BED specification is at http://www.genome.ucsc.edu/FAQ/FAQformat.html
   The VCF specification is at https://github.com/samtools/hts-specs
 */

void *bed_read(const char *fn)
{
    reghash_t *h = kh_init(reg);
    gzFile fp;
    kstream_t *ks = NULL;
    int dret;
    unsigned int line = 0;
    kstring_t str = { 0, 0, NULL };

    if (NULL == h) return NULL;
    // read the list
    fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) return 0;
    ks = ks_init(fp);
    if (NULL == ks) goto fail;  // In case ks_init ever gets error checking...
    while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) > 0) { // read a line
        char *ref = str.s, *ref_end;
        unsigned int beg = 0, end = 0;
        int num = 0;
        khint_t k;
        bed_reglist_t *p;

        line++;
        while (*ref && isspace(*ref)) ref++;
        if ('\0' == *ref) continue;  // Skip blank lines
        if ('#'  == *ref) continue;  // Skip BED file comments
        ref_end = ref;   // look for the end of the reference name
        while (*ref_end && !isspace(*ref_end)) ref_end++;
        if ('\0' != *ref_end) {
            *ref_end = '\0';  // terminate ref and look for start, end
            num = sscanf(ref_end + 1, "%u %u", &beg, &end);
        }
        if (1 == num) {  // VCF-style format
            end = beg--; // Counts from 1 instead of 0 for BED files
        }
        if (num < 1 || end < beg) {
            // These two are special lines that can occur in BED files.
            // Check for them here instead of earlier in case someone really
            // has called their reference "browser" or "track".
            if (0 == strcmp(ref, "browser")) continue;
            if (0 == strcmp(ref, "track")) continue;
            fprintf(stderr, "[bed_read] Parse error reading %s at line %u\n",
                    fn, line);
            goto fail_no_msg;
        }

        // Put reg in the hash table if not already there
        k = kh_get(reg, h, ref);
        if (k == kh_end(h)) { // absent from the hash table
            int ret;
            char *s = strdup(ref);
            if (NULL == s) goto fail;
            k = kh_put(reg, h, s, &ret);
            if (-1 == ret) {
                free(s);
                goto fail;
            }
            memset(&kh_val(h, k), 0, sizeof(bed_reglist_t));
        }
        p = &kh_val(h, k);

        // Add begin,end to the list
        if (p->n == p->m) {
            p->m = p->m? p->m<<1 : 4;
            p->a = realloc(p->a, p->m * 8);
            if (NULL == p->a) goto fail;
        }
        p->a[p->n++] = (uint64_t)beg<<32 | end;
    }
    // FIXME: Need to check for errors in ks_getuntil.  At the moment it
    // doesn't look like it can return one.  Possibly use gzgets instead?

    ks_destroy(ks);
    gzclose(fp);
    free(str.s);
    bed_index(h);
    return h;
 fail:
    fprintf(stderr, "[bed_read] Error reading %s : %s\n", fn, strerror(errno));
 fail_no_msg:
    if (ks) ks_destroy(ks);
    if (fp) gzclose(fp);
    free(str.s);
    bed_destroy(h);
    return NULL;
}

void bed_destroy(void *_h)
{
    reghash_t *h = (reghash_t*)_h;
    khint_t k;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free(kh_val(h, k).a);
            free(kh_val(h, k).idx);
            free((char*)kh_key(h, k));
        }
    }
    kh_destroy(reg, h);
}
