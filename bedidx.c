/*  bedidx.c -- BED file indexing.

    Copyright (C) 2011 Broad Institute.
    Copyright (C) 2014, 2017-2019 Genome Research Ltd.

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
#include "bedidx.h"

#include "htslib/ksort.h"

#include "htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 8192)

static inline int lt_pair_pos(hts_pair_pos_t a, hts_pair_pos_t b) {
    if (a.beg == b.beg) return a.end < b.end;
    return a.beg < b.beg;
}
KSORT_INIT_STATIC(hts_pair_pos_t, hts_pair_pos_t, lt_pair_pos)

/*! @typedef
 * @abstract bed_reglist_t - value type of the BED hash table
 * This structure encodes the list of intervals (ranges) for the regions provided via BED file or
 * command line arguments.
 * @field *a           pointer to the array of intervals.
 * @field n            actual number of elements contained by a
 * @field m            number of allocated elements to a (n <= m)
 * @field *idx         index array for computing the minimum offset
 */
typedef struct {
    int n, m;
    hts_pair_pos_t *a;
    int *idx;
    int filter;
} bed_reglist_t;

#include "htslib/khash.h"
KHASH_MAP_INIT_STR(reg, bed_reglist_t)

typedef kh_reg_t reghash_t;

#if 0
// Debug function
static void bed_print(void *reg_hash) {
    reghash_t *h = (reghash_t *)reg_hash;
    bed_reglist_t *p;
    khint_t k;
    int i;
    const char *reg;

    if (!h) {
        printf("Hash table is empty!\n");
        return;
    }
    for (k = kh_begin(h); k < kh_end(h); k++) {
        if (kh_exist(h,k)) {
            reg = kh_key(h,k);
            printf("Region: '%s'\n", reg);
            if ((p = &kh_val(h,k)) != NULL && p->n > 0) {
                printf("Filter: %d\n", p->filter);
                for (i=0; i<p->n; i++) {
                    printf("\tinterval[%d]: %"PRIhts_pos"-%"PRIhts_pos"\n",
                           i,p->a[i].beg,p->a[i].end);
                }
            } else {
                printf("Region '%s' has no intervals!\n", reg);
            }
        }
    }
}
#endif

static int *bed_index_core(int n, hts_pair_pos_t *a)
{
    int i, j, l, *idx, *new_idx;
    l = 0; idx = 0;
    for (i = 0; i < n; ++i) {
        hts_pos_t beg, end;
        beg = a[i].beg >> LIDX_SHIFT; end = a[i].end >> LIDX_SHIFT;
        if (l < end + 1) {
            int old_l = l;
            l = end + 1;
            kroundup32(l);
            new_idx = realloc(idx, l * sizeof(*idx));
            if (!new_idx) {
                free(idx);
                return NULL;
            }
            idx = new_idx;

            for (j = old_l; j < l; ++j)
                idx[j] = -1;
        }

        for (j = beg; j < end+1; ++j)
            if (idx[j] < 0)
                idx[j] = i;
    }
    return idx;
}

static void bed_index(void *_h)
{
    reghash_t *h = (reghash_t*)_h;
    khint_t k;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            bed_reglist_t *p = &kh_val(h, k);
            if (p->idx) free(p->idx);
            ks_introsort(hts_pair_pos_t, p->n, p->a);
            p->idx = bed_index_core(p->n, p->a);
        }
    }
}

static int bed_minoff(const bed_reglist_t *p, hts_pos_t beg, hts_pos_t end) {
    int i, min_off=0;

    if (p && p->idx) {
        min_off = (beg>>LIDX_SHIFT >= p->n)? p->idx[p->n-1] : p->idx[beg>>LIDX_SHIFT];
        if (min_off < 0) { // TODO: this block can be improved, but speed should not matter too much here
            hts_pos_t n = beg>>LIDX_SHIFT;
            if (n > p->n)
                n = p->n;
            for (i = n - 1; i >= 0; --i)
                if (p->idx[i] >= 0)
                    break;
            min_off = i >= 0? p->idx[i] : 0;
        }
    }

    return min_off;
}

static int bed_overlap_core(const bed_reglist_t *p, hts_pos_t beg, hts_pos_t end)
{
    int i, min_off;
    if (p->n == 0) return 0;
    min_off = bed_minoff(p, beg, end);

    for (i = min_off; i < p->n; ++i) {
        if (p->a[i].beg >= end) break; // out of range; no need to proceed
        if (p->a[i].end > beg && p->a[i].beg < end)
            return 1; // find the overlap; return
    }
    return 0;
}

int bed_overlap(const void *_h, const char *chr, hts_pos_t beg, hts_pos_t end)
{
    const reghash_t *h = (const reghash_t*)_h;
    khint_t k;
    if (!h) return 0;
    k = kh_get(reg, h, chr);
    if (k == kh_end(h)) return 0;
    return bed_overlap_core(&kh_val(h, k), beg, end);
}

/** @brief Trim a sorted interval list, inside a region hash table,
 *   by removing completely contained intervals and merging adjacent or
 *   overlapping intervals.
 *  @param reg_hash    the region hash table with interval lists as values
 */

void bed_unify(void *reg_hash) {

    int i, j, new_n;
    reghash_t *h;
    bed_reglist_t *p;

    if (!reg_hash)
        return;

    h = (reghash_t *)reg_hash;

    for (i = kh_begin(h); i < kh_end(h); i++) {
        if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || !(p->n))
            continue;

        for (new_n = 0, j = 1; j < p->n; j++) {
            if (p->a[new_n].end < p->a[j].beg) {
                p->a[++new_n] = p->a[j];
            } else {
                if (p->a[new_n].end < p->a[j].end)
                    p->a[new_n].end = p->a[j].end;
            }
        }

        p->n = ++new_n;
    }
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
    unsigned int line = 0, save_errno;
    kstring_t str = { 0, 0, NULL };

    if (NULL == h) return NULL;
    // read the list
    fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) return 0;
    ks = ks_init(fp);
    if (NULL == ks) goto fail;  // In case ks_init ever gets error checking...
    int ks_len;
    while ((ks_len = ks_getuntil(ks, KS_SEP_LINE, &str, &dret)) >= 0) { // read a line
        char *ref = str.s, *ref_end;
        uint64_t beg = 0, end = 0;
        int num = 0;
        khint_t k;
        bed_reglist_t *p;

        if (ks_len == 0)
            continue; // skip blank lines

        line++;
        while (*ref && isspace(*ref)) ref++;
        if ('\0' == *ref) continue;  // Skip blank lines
        if ('#'  == *ref) continue;  // Skip BED file comments
        ref_end = ref;   // look for the end of the reference name
        while (*ref_end && !isspace(*ref_end)) ref_end++;
        if ('\0' != *ref_end) {
            *ref_end = '\0';  // terminate ref and look for start, end
            num = sscanf(ref_end + 1, "%"SCNu64" %"SCNu64, &beg, &end);
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
            if (num < 1) {
                fprintf(stderr,
                        "[bed_read] Parse error reading \"%s\" at line %u\n",
                        fn, line);
            } else {
                fprintf(stderr,
                        "[bed_read] Parse error reading \"%s\" at line %u : "
                        "end (%"PRIu64") must not be less "
                        "than start (%"PRIu64")\n",
                        fn, line, end, beg);
            }
            errno = 0; // Prevent caller from printing misleading error messages
            goto fail;
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
            p->m = p->m ? p->m<<1 : 4;
            hts_pair_pos_t *new_a = realloc(p->a, p->m * sizeof(p->a[0]));
            if (NULL == new_a) goto fail;
            p->a = new_a;
        }
        p->a[p->n].beg = beg;
        p->a[p->n++].end = end;
    }
    // FIXME: Need to check for errors in ks_getuntil.  At the moment it
    // doesn't look like it can return one.  Possibly use gzgets instead?

    if (gzclose(fp) != Z_OK) {
        fp = NULL;
        goto fail;
    }
    ks_destroy(ks);
    free(str.s);
    bed_index(h);
    //bed_unify(h);
    return h;
 fail:
    save_errno = errno;
    if (ks) ks_destroy(ks);
    if (fp) gzclose(fp);
    free(str.s);
    bed_destroy(h);
    errno = save_errno;
    return NULL;
}

void bed_destroy(void *_h)
{
    reghash_t *h;
    khint_t k;

    if (!_h)
        return;

    h = (reghash_t*)_h;
    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free(kh_val(h, k).a);
            free(kh_val(h, k).idx);
            free((char*)kh_key(h, k));
        }
    }
    kh_destroy(reg, h);
}

static void *bed_insert(void *reg_hash, char *reg, hts_pos_t beg, hts_pos_t end) {

    reghash_t *h;
    khint_t k;
    bed_reglist_t *p;

    if (!reg_hash)
        return NULL;

    h = (reghash_t *)reg_hash;

    // Put reg in the hash table if not already there
    k = kh_get(reg, h, reg); //looks strange, but only the second reg is the actual region name.
    if (k == kh_end(h)) { // absent from the hash table
        int ret;
        char *s = strdup(reg);
        if (NULL == s) goto fail;
        k = kh_put(reg, h, s, &ret);
        if (-1 == ret) {
            free(s);
            goto fail;
        }
        memset(&kh_val(h, k), 0, sizeof(bed_reglist_t));
    }
    p = &kh_val(h, k);

    // Add beg and end to the list
    if (p->n == p->m) {
        p->m = p->m ? p->m<<1 : 4;
        hts_pair_pos_t *new_a = realloc(p->a, p->m * sizeof(p->a[0]));
        if (NULL == new_a) goto fail;
        p->a = new_a;
    }
    p->a[p->n].beg = beg;
    p->a[p->n++].end = end;

fail:
    return h;
}

/* @brief Filter a region hash table (coming from the BED file) by another
 *  region hash table (coming from CLI), so that only intervals contained in
 *  both hash tables are kept.
 * @param reg_hash    the target region hash table
 * @param tmp_hash    the filter region hash table
 * @return            pointer to the filtered hash table
 */

static void *bed_filter(void *reg_hash, void *tmp_hash) {

    reghash_t *h;
    reghash_t *t;
    bed_reglist_t *p, *q;
    khint_t l, k;
    hts_pair_pos_t *new_a;
    int i, j, new_n, min_off;
    const char *reg;
    hts_pos_t beg, end;

    h = (reghash_t *)reg_hash;
    t = (reghash_t *)tmp_hash;
    if (!h)
        return NULL;
    if (!t)
        return h;

    for (l = kh_begin(t); l < kh_end(t); l++) {
        if (!kh_exist(t,l) || !(q = &kh_val(t,l)) || !(q->n))
            continue;

        reg = kh_key(t,l);
        k = kh_get(reg, h, reg); //looks strange, but only the second reg is a proper argument.
        if (k == kh_end(h) || !(p = &kh_val(h, k)) || !(p->n))
            continue;

        new_a = calloc(q->n + p->n, sizeof(new_a[0]));
        if (!new_a)
            return NULL;
        new_n = 0;

        for (i = 0; i < q->n; i++) {
            beg = q->a[i].beg;
            end = q->a[i].end;

            min_off = bed_minoff(p, beg, end);
            for (j = min_off; j < p->n; ++j) {
                if (p->a[j].beg >= end) break; // out of range; no need to proceed
                if (p->a[j].end > beg && p->a[j].beg < end) {
                    new_a[new_n].beg = MAX(p->a[j].beg, beg);
                    new_a[new_n++].end = MIN(p->a[j].end, end);
                }
            }
        }

        if (new_n > 0) {
            free(p->a);
            p->a = new_a;
            p->n = new_n;
            p->m = new_n;
            p->filter = FILTERED;
        } else {
            free(new_a);
            p->filter = ALL;
        }
    }

    return h;
}

void *bed_hash_regions(void *reg_hash, char **regs, int first, int last, int *op) {

    reghash_t *h = (reghash_t *)reg_hash;
    reghash_t *t = NULL;

    int i;
    char reg[1024];
    const char *q;
    int beg, end;

    if (h) {
        t = kh_init(reg);
        if (!t) {
            fprintf(stderr, "Error when creating the temporary region hash table!\n");
            return NULL;
        }
    } else {
        h = kh_init(reg);
        if (!h) {
            fprintf(stderr, "Error when creating the region hash table!\n");
            return NULL;
        }
        *op = 1;
    }

    for (i=first; i<last; i++) {

        // Note, ideally we would call sam_parse_region here, but it's complicated by not
        // having the sam header known and the likelihood of the bed file containing data for other
        // references too which we currently just ignore.
        //
        // TO DO...
        q = hts_parse_reg(regs[i], &beg, &end);
        if (q) {
            if ((int)(q - regs[i] + 1) > 1024) {
                fprintf(stderr, "Region name '%s' is too long (bigger than %d).\n", regs[i], 1024);
                continue;
            }
            strncpy(reg, regs[i], q - regs[i]);
            reg[q - regs[i]] = 0;
        } else {
            // not parsable as a region, but possibly a sequence named "foo:a"
            if (strlen(regs[i]) + 1 > 1024) {
                fprintf(stderr, "Region name '%s' is too long (bigger than %d).\n", regs[i], 1024);
                continue;
            }
            strcpy(reg, regs[i]);
            beg = 0; end = INT_MAX;
        }

        //if op==1 insert reg to the bed hash table
        if (*op && !(bed_insert(h, reg, beg, end))) {
            fprintf(stderr, "Error when inserting region='%s' in the bed hash table at address=%p!\n", regs[i], h);
        }
        //if op==0, first insert the regions in the temporary hash table,
        //then filter the bed hash table using it
        if (!(*op) && !(bed_insert(t, reg, beg, end))) {
            fprintf(stderr, "Error when inserting region='%s' in the temporary hash table at address=%p!\n", regs[i], t);
        }
    }

    if (!(*op)) {
        bed_index(t);
        bed_unify(t);
        h = bed_filter(h, t);
        bed_destroy(t);
    }

    if (h) {
        bed_index(h);
        bed_unify(h);
    }

    return h;
}

const char* bed_get(void *reg_hash, int i, int filter) {

    reghash_t *h;
    bed_reglist_t *p;

    if (!reg_hash)
        return NULL;

    h = (reghash_t *)reg_hash;
    if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || (p->filter < filter))
        return NULL;

    return kh_key(h, i);
}

hts_reglist_t *bed_reglist(void *reg_hash, int filter, int *n_reg) {

    reghash_t *h;
    bed_reglist_t *p;
    khint_t i;
    hts_reglist_t *reglist = NULL;
    int count = 0;
    int j;

    if (!reg_hash)
        return NULL;

    h = (reghash_t *)reg_hash;

    for (i = kh_begin(h); i < kh_end(h); i++) {
        if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || (p->filter < filter))
            continue;
        count++;
    }
    if (!count)
        return NULL;

    reglist = (hts_reglist_t *)calloc(count, sizeof(hts_reglist_t));
    if (!reglist)
        return NULL;

    *n_reg = count;
    count = 0;

    for (i = kh_begin(h); i < kh_end(h) && count < *n_reg; i++) {
        if (!kh_exist(h,i) || !(p = &kh_val(h,i)) || (p->filter < filter))
            continue;

        reglist[count].reg = kh_key(h,i);
        reglist[count].intervals = (hts_pair32_t *)calloc(p->n, sizeof(hts_pair32_t));
        if(!(reglist[count].intervals)) {
            hts_reglist_free(reglist, count);
            return NULL;
        }
        reglist[count].count = p->n;
        reglist[count].max_end = 0;

        for (j = 0; j < p->n; j++) {
            reglist[count].intervals[j].beg = p->a[j].beg;
            reglist[count].intervals[j].end = p->a[j].end;

            if (reglist[count].intervals[j].end > reglist[count].max_end)
                reglist[count].max_end = reglist[count].intervals[j].end;
        }
        count++;
    }

    return reglist;
}
