/*  sam_view.c -- SAM<->BAM<->CRAM conversion.

    Copyright (C) 2009-2023 Genome Research Ltd.
    Portions copyright (C) 2009, 2011, 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
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
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "htslib/hts_expr.h"
#include "samtools.h"
#include "sam_opts.h"
#include "bam.h" // for bam_get_library and bam_remove_B
#include "bedidx.h"
#include "sam_utils.h"

KHASH_SET_INIT_STR(str)
typedef khash_t(str) *strhash_t;

// This structure contains the settings for a samview run
typedef struct samview_settings {
    strhash_t rghash;
    strhash_t rnhash;
    strhash_t tvhash;
    int min_mapQ;

    // Described here in the same terms as the usage statement.
    // The code however always negates to "reject if"         keep if:
    int flag_on;     // keep   if (FLAG & N) == N             (all on)
    int flag_off;    // keep   if (FLAG & N) == 0             (all off)
    int flag_anyon;  // keep   if (FLAG & N) != 0             (any on)
    int flag_alloff; // reject if (FLAG & N) == N             (any off)

    int min_qlen;
    int remove_B;
    uint32_t subsam_seed;
    double subsam_frac;
    char* library;
    void* bed;
    size_t remove_aux_len;
    char** remove_aux;
    int multi_region;
    char* tag;
    hts_filter_t *filter;
    int remove_flag;
    int add_flag;
    int unmap;
    auxhash_t remove_tag;
    auxhash_t keep_tag;

    hts_idx_t *hts_idx;
    sam_hdr_t *header;
    samFile *in, *out, *un_out;
    int64_t count;
    int is_count;
    char *fn_in, *fn_idx_in, *fn_out, *fn_fai, *fn_un_out, *fn_out_idx, *fn_un_out_idx;
    int fetch_pairs, nreglist;
    hts_reglist_t *reglist;
    int sanitize;
    int count_rf; // CRAM_OPT_REQUIRED_FIELDS for view -c
} samview_settings_t;

// Copied from htslib/sam.c.
// TODO: we need a proper interface to find the length of an aux tag,
// or at the very make exportable versions of these in htslib.
static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

// Copied from htslib/sam.c.
static inline uint8_t *skip_aux(uint8_t *s, uint8_t *end)
{
    int size;
    uint32_t n;
    if (s >= end) return end;
    size = aux_type2size(*s); ++s; // skip type
    switch (size) {
    case 'Z':
    case 'H':
        while (s < end && *s) ++s;
        return s < end ? s + 1 : end;
    case 'B':
        if (end - s < 5) return NULL;
        size = aux_type2size(*s); ++s;
        n = le_to_u32(s);
        s += 4;
        if (size == 0 || end - s < size * n) return NULL;
        return s + size * n;
    case 0:
        return NULL;
    default:
        if (end - s < size) return NULL;
        return s + size;
    }
}

// Returns 0 to indicate read should be output 1 otherwise
static int process_aln(const sam_hdr_t *h, bam1_t *b, samview_settings_t* settings)
{
    if (settings->filter && sam_passes_filter(h, b, settings->filter) < 1)
        return 1;

    if (settings->remove_B) bam_remove_B(b);
    if (settings->min_qlen > 0) {
        int k, qlen = 0;
        uint32_t *cigar = bam_get_cigar(b);
        for (k = 0; k < b->core.n_cigar; ++k)
            if ((bam_cigar_type(bam_cigar_op(cigar[k]))&1) || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
                qlen += bam_cigar_oplen(cigar[k]);
        if (qlen < settings->min_qlen) return 1;
    }
    if (b->core.qual < settings->min_mapQ || ((b->core.flag & settings->flag_on) != settings->flag_on) || (b->core.flag & settings->flag_off))
        return 1;
    if (settings->flag_alloff && ((b->core.flag & settings->flag_alloff) == settings->flag_alloff))
        return 1;
    if (settings->flag_anyon && ((b->core.flag & settings->flag_anyon) == 0))
        return 1;
    if (!settings->multi_region && settings->bed && (b->core.tid < 0 || !bed_overlap(settings->bed, sam_hdr_tid2name(h, b->core.tid), b->core.pos, bam_endpos(b))))
        return 1;
    if (settings->subsam_frac > 0.) {
        uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ settings->subsam_seed);
        if ((double)(k&0xffffff) / 0x1000000 >= settings->subsam_frac) return 1;
    }
    if (settings->rghash) {
        uint8_t *s = bam_aux_get(b, "RG");
        if (s) {
            khint_t k = kh_get(str, settings->rghash, (char*)(s + 1));
            if (k == kh_end(settings->rghash)) return 1;
        }
    }
    if (settings->tag) {
        uint8_t *s = bam_aux_get(b, settings->tag);
        if (s) {
            if (settings->tvhash) {
                char t[32], *val;
                if (*s == 'i' || *s == 'I' || *s == 's' || *s == 'S' || *s == 'c' || *s == 'C') {
                    int ret = snprintf(t, 32, "%"PRId64, bam_aux2i(s));
                    if (ret > 0) val = t;
                    else return 1;
                } else if (*s == 'A') {
                    t[0] = *(s+1);
                    t[1] = 0;
                    val = t;
                } else {
                    val = (char *)(s+1);
                }
                khint_t k = kh_get(str, settings->tvhash, val);
                if (k == kh_end(settings->tvhash)) return 1;
            }
        } else {
            return 1;
        }
    }
    if (settings->rnhash) {
        const char* rn = bam_get_qname(b);
        if (!rn || kh_get(str, settings->rnhash, rn) == kh_end(settings->rnhash)) {
            return 1;
        }
    }
    if (settings->library) {
        const char *p = bam_get_library((sam_hdr_t*)h, b);
        if (!p || strcmp(p, settings->library) != 0) return 1;
    }
    return 0;
}

static int adjust_tags(const sam_hdr_t *h, bam1_t *b,
                       samview_settings_t* settings) {
    if (settings->keep_tag) {
        uint8_t *s_from, *s_to, *end = b->data + b->l_data;
        auxhash_t h = settings->keep_tag;

        s_from = s_to = bam_get_aux(b);
        while (s_from < end) {
            int x = (int)s_from[0]<<8 | s_from[1];
            uint8_t *s = skip_aux(s_from+2, end);
            if (s == NULL) {
                print_error("view", "malformed aux data for record \"%s\"",
                            bam_get_qname(b));
                return -1;
            }

            if (kh_get(aux_exists, h, x) != kh_end(h) ) {
                if (s_to != s_from) memmove(s_to, s_from, s - s_from);
                s_to += s - s_from;
            }
            s_from = s;
        }
        b->l_data = s_to - b->data;

    } else if (settings->remove_tag) {
        uint8_t *s_from, *s_to, *end = b->data + b->l_data;
        auxhash_t h = settings->remove_tag;

        s_from = s_to = bam_get_aux(b);
        while (s_from < end) {
            int x = (int)s_from[0]<<8 | s_from[1];
            uint8_t *s = skip_aux(s_from+2, end);
            if (s == NULL) {
                print_error("view", "malformed aux data for record \"%s\"",
                            bam_get_qname(b));
                return -1;
            }

            if (kh_get(aux_exists, h, x) == kh_end(h) ) {
                if (s_to != s_from) memmove(s_to, s_from, s - s_from);
                s_to += s - s_from;
            }
            s_from = s;
        }
        b->l_data = s_to - b->data;
    }

    return 0;
}

static int usage(FILE *fp, int exit_status, int is_long_help);

static int populate_lookup_from_file(const char *subcmd, strhash_t lookup, char *fn)
{
    FILE *fp;
    char buf[1024];
    int ret = 0;
    fp = fopen(fn, "r");
    if (fp == NULL) {
        print_error_errno(subcmd, "failed to open \"%s\" for reading", fn);
        return -1;
    }

    while (ret != -1 && !feof(fp) && fscanf(fp, "%1023s", buf) > 0) {
        char *d = strdup(buf);
        if (d != NULL) {
            kh_put(str, lookup, d, &ret);
            if (ret == 0) free(d); /* Duplicate */
        } else {
            ret = -1;
        }
    }
    if (ferror(fp)) ret = -1;
    if (ret == -1) {
        print_error_errno(subcmd, "failed to read \"%s\"", fn);
    }
    fclose(fp);
    return (ret != -1) ? 0 : -1;
}

static int add_read_group_single(const char *subcmd, samview_settings_t *settings, char *name)
{
    char *d = strdup(name);
    int ret = 0;

    if (d == NULL) goto err;

    if (settings->rghash == NULL) {
        settings->rghash = kh_init(str);
        if (settings->rghash == NULL) goto err;
    }

    kh_put(str, settings->rghash, d, &ret);
    if (ret == -1) goto err;
    if (ret ==  0) free(d); /* Duplicate */
    return 0;

 err:
    print_error(subcmd, "Couldn't add \"%s\" to read group list: memory exhausted?", name);
    free(d);
    return -1;
}

static int add_read_names_file(const char *subcmd, samview_settings_t *settings, char *fn)
{
    if (settings->rnhash == NULL) {
        settings->rnhash = kh_init(str);
        if (settings->rnhash == NULL) {
            perror(NULL);
            return -1;
        }
    }
    return populate_lookup_from_file(subcmd, settings->rnhash, fn);
}

static int add_read_groups_file(const char *subcmd, samview_settings_t *settings, char *fn)
{
    if (settings->rghash == NULL) {
        settings->rghash = kh_init(str);
        if (settings->rghash == NULL) {
            perror(NULL);
            return -1;
        }
    }
    return populate_lookup_from_file(subcmd, settings->rghash, fn);
}

static int add_tag_value_single(const char *subcmd, samview_settings_t *settings, char *name)
{
    char *d = strdup(name);
    int ret = 0;

    if (d == NULL) goto err;

    if (settings->tvhash == NULL) {
        settings->tvhash = kh_init(str);
        if (settings->tvhash == NULL) goto err;
    }

    kh_put(str, settings->tvhash, d, &ret);
    if (ret == -1) goto err;
    if (ret ==  0) free(d); /* Duplicate */
    return 0;

 err:
    print_error(subcmd, "Couldn't add \"%s\" to tag values list: memory exhausted?", name);
    free(d);
    return -1;
}

static int add_tag_values_file(const char *subcmd, samview_settings_t *settings, char *fn)
{
    if (settings->tvhash == NULL) {
        settings->tvhash = kh_init(str);
        if (settings->tvhash == NULL) {
            perror(NULL);
            return -1;
        }
    }
    return populate_lookup_from_file(subcmd, settings->tvhash, fn);
}

static inline int check_sam_write1(samFile *fp, const sam_hdr_t *h, const bam1_t *b, const char *fname, int *retp)
{
    int r = sam_write1(fp, h, b);
    if (r >= 0) return r;

    if (fname) print_error_errno("view", "writing to \"%s\" failed", fname);
    else print_error_errno("view", "writing to standard output failed");

    *retp = EXIT_FAILURE;
    return r;
}

static inline void change_flag(bam1_t *b, samview_settings_t *settings)
{
    if (settings->add_flag)
        b->core.flag |= settings->add_flag;

    if (settings->remove_flag)
        b->core.flag &= ~settings->remove_flag;
}

static int cmp_reglist_intervals(const void *aptr, const void *bptr)
{
    hts_pair_pos_t *a = (hts_pair_pos_t*)aptr;
    hts_pair_pos_t *b = (hts_pair_pos_t*)bptr;
    if ( a->beg < b->beg ) return -1;
    if ( a->beg > b->beg ) return 1;
    if ( a->end < b->end ) return -1;
    if ( a->end > b->end ) return 1;
    return 0;
}
static int cmp_reglist_tids(const void *aptr, const void *bptr)
{
    hts_reglist_t *a = (hts_reglist_t*)aptr;
    hts_reglist_t *b = (hts_reglist_t*)bptr;
    if ( b->tid==HTS_IDX_NOCOOR || a->tid < b->tid ) return -1;
    if ( a->tid==HTS_IDX_NOCOOR || a->tid > b->tid ) return 1;
    return 0;
}

static hts_reglist_t *_reglist_dup(sam_hdr_t *hdr, hts_reglist_t *src, int nsrc)
{
    int i,j;
    hts_reglist_t *dst = (hts_reglist_t*)calloc(nsrc,sizeof(hts_reglist_t));
    if ( !dst ) {
        print_error_errno("view", "[%s:%d] could not allocate region list"
                          ,__FILE__ ,__LINE__);
        return NULL;
    }
    for (i=0; i<nsrc; i++)
    {
        // Assume tid is not set correctly, reg is informative but may not point to a long-lived memory
        dst[i].tid = sam_hdr_name2tid(hdr,src[i].reg);
        dst[i].min_beg = src[i].min_beg;
        dst[i].max_end = src[i].max_end;
        dst[i].count = src[i].count;
        dst[i].intervals = (hts_pair_pos_t*)malloc(sizeof(hts_pair_pos_t)*dst[i].count);
        if ( !dst[i].intervals ) {
            print_error_errno("view", "[%s:%d] could not allocate region list",
                              __FILE__, __LINE__);
            goto fail;
        }
        for (j=0; j<dst[i].count; j++)
            dst[i].intervals[j] = src[i].intervals[j];
    }
    qsort(dst,nsrc,sizeof(*dst),cmp_reglist_tids);
    return dst;

 fail:
    for (j = 0; j < i; j++)
        free(dst[j].intervals);
    free(dst);
    return NULL;
}
static inline int _reglist_find_tid(hts_reglist_t *reg, int nreg, int tid) // binary search
{
    int i = -1, imin = 0, imax = nreg - 1;
    while ( imin <= imax )
    {
        i = (imin+imax)/2;
        if ( tid==HTS_IDX_NOCOOR || reg[i].tid < tid ) imin = i + 1;
        else if ( reg[i].tid==HTS_IDX_NOCOOR || reg[i].tid > tid ) imax = i - 1;
        else break;
    }
    if ( i<0 || reg[i].tid < tid ) i++;    // not found, i will be the index of the inserted element
    return i;
}
static int _reglist_push(hts_reglist_t **_reg, int *_nreg, int tid, hts_pos_t beg, hts_pos_t end)
{
    hts_reglist_t *reg = *_reg;
    int nreg = *_nreg;
    int i = _reglist_find_tid(reg,nreg,tid);
    if ( i>=nreg || reg[i].tid!=tid ) {
        nreg++;
        reg = (hts_reglist_t*)realloc(reg,sizeof(hts_reglist_t)*nreg);
        if ( !reg ) {
            print_error_errno("view", "[%s:%d] could not extend region list",
                              __FILE__, __LINE__);
            return -1;
        }
        if ( i+1 < nreg )
            memmove(reg + i + 1, reg + i, sizeof(hts_reglist_t)*(nreg - i - 1));
        reg[i].reg = NULL;
        reg[i].tid = tid;
        reg[i].min_beg = beg;
        reg[i].max_end = end;
        reg[i].intervals = NULL;
        reg[i].count = 0;
    }
    *_reg = reg;
    *_nreg = nreg;
    if ( reg[i].count > 0
         && reg[i].intervals[reg[i].count - 1].beg==beg
         && reg[i].intervals[reg[i].count - 1].end==end ) {
        return 0;
    }
    hts_pair_pos_t *new_intervals = realloc(reg[i].intervals, sizeof(hts_pair_pos_t)*(reg[i].count + 1));
    if (!new_intervals) {
        print_error_errno("view", "[%s:%d] could not extend region list",
                          __FILE__, __LINE__);
        return -1;
    }
    reg[i].intervals = new_intervals;
    reg[i].intervals[reg[i].count].beg = beg;
    reg[i].intervals[reg[i].count].end = end;
    reg[i].count++;
    return 0;
}

static void _reglist_merge(hts_reglist_t *reg, int nreg)
{
    int i,j;
    for (i=0; i<nreg; i++)
    {
        qsort(reg[i].intervals,reg[i].count,sizeof(*reg[i].intervals),cmp_reglist_intervals);
        int k = 1;
        for (j=1; j<reg[i].count; j++)
        {
            if ( reg[i].intervals[k-1].end < reg[i].intervals[j].beg )
            {
                if ( k < j ) reg[i].intervals[k] = reg[i].intervals[j];
                k++;
                continue;
            }
            if ( reg[i].intervals[k-1].end < reg[i].intervals[j].end ) reg[i].intervals[k-1].end = reg[i].intervals[j].end;
        }
        reg[i].count = k;
        reg[i].max_end = reg[i].intervals[k-1].end;
    }
}
hts_itr_multi_t *multi_region_init(samview_settings_t *conf, char **regs, int nregs)
{
    hts_itr_multi_t *iter = NULL;
    int filter_state = ALL;
    if ( nregs ) {
        int filter_op = 0;
        conf->bed = bed_hash_regions(conf->bed, regs, 0, nregs, &filter_op); // insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
        if ( !filter_op )
            filter_state = FILTERED;
    }
    else
        bed_unify(conf->bed);
    if ( !conf->bed) { // index is unavailable or no regions have been specified
        print_error("view", "No regions or BED file have been provided. Aborting.");
        return NULL;
    }

    int regcount = 0;
    hts_reglist_t *reglist = bed_reglist(conf->bed, filter_state, &regcount);
    if (!reglist) {
        print_error("view", "Region list is empty or could not be created. Aborting.");
        return NULL;
    }

    if ( conf->fetch_pairs ) {
        conf->reglist  = _reglist_dup(conf->header,reglist,regcount);
        if (!conf->reglist)
            return NULL;
        conf->nreglist = regcount;
    }

    iter = sam_itr_regions(conf->hts_idx, conf->header, reglist, regcount);
    if ( !iter ) {
        print_error("view", "Iterator could not be created. Aborting.");
        return NULL;
    }
    return iter;
}

KHASH_SET_INIT_STR(names)

static int fetch_pairs_collect_mates(samview_settings_t *conf, hts_itr_multi_t *iter)
{
    khint_t k;
    int nunmap = 0, r = 0, nmates = 0, write_error = 0, retval = EXIT_FAILURE;
    kh_names_t *mate_names = kh_init(names);
    bam1_t *rec = bam_init1();

    if (!mate_names) {
        print_error_errno("view", "could not allocate mate names table");
        goto out;
    }
    if (!rec) {
        print_error_errno("view", "could not allocate bam record");
        goto out;
    }

    while ((r =sam_itr_multi_next(conf->in, iter, rec))>=0) {
        if ( (rec->core.flag & BAM_FPAIRED) == 0 ) continue;
        if ( rec->core.mtid>=0 && bed_overlap(conf->bed, sam_hdr_tid2name(conf->header,rec->core.mtid), rec->core.mpos, rec->core.mpos) ) continue;
        if ( process_aln(conf->header, rec, conf) ) continue;

        nmates++;

        k = kh_get(names,mate_names,bam_get_qname(rec));
        if ( k == kh_end(mate_names) ) {
            int ret = 0;
            char *name_copy = strdup(bam_get_qname(rec));
            if (!name_copy) {
                print_error_errno("view", "[%s:%d] could not store sample name, %d elements", __FILE__,__LINE__,nmates);
                goto out;
            }
            kh_put(names, mate_names, name_copy, &ret);
            if ( ret<0 ) {
                print_error_errno("view", "[%s:%d] could not store sample name, %d elements",__FILE__,__LINE__,nmates);
                free(name_copy);
                goto out;
            }
        }

        if ( rec->core.mtid < 0 || (rec->core.flag & BAM_FMUNMAP) ) nunmap = 1;
        if ( rec->core.mtid >= 0 ) {
            if (_reglist_push(&conf->reglist, &conf->nreglist, rec->core.mtid, rec->core.mpos,rec->core.mpos+1) != 0)
                goto out;
        }
    }

    if (r < -1) {
        print_error_errno("view", "error reading file \"%s\"", conf->fn_in);
        goto out;
    }

    _reglist_merge(conf->reglist, conf->nreglist);
    if ( nunmap ) {
        if (_reglist_push(&conf->reglist,&conf->nreglist,HTS_IDX_NOCOOR,0,HTS_POS_MAX) != 0)
            goto out;
    }
    hts_itr_multi_destroy(iter);
    iter = sam_itr_regions(conf->hts_idx, conf->header, conf->reglist, conf->nreglist);
    if ( !iter ) {
        print_error_errno("view", "[%s:%d] iterator could not be created",__FILE__,__LINE__);
        goto out;
    }
    while ((r = sam_itr_multi_next(conf->in, iter, rec))>=0) {
        int drop = 1;
        if (rec->core.tid >=0 &&
            bed_overlap(conf->bed, sam_hdr_tid2name(conf->header,rec->core.tid), rec->core.pos, bam_endpos(rec))) drop = 0;
        if ( drop ) {
             k = kh_get(names,mate_names,bam_get_qname(rec));
             if ( k != kh_end(mate_names) ) drop = 0;
        }
        if (!drop && process_aln(conf->header, rec, conf) == 0) {
            if (adjust_tags(conf->header, rec, conf) != 0)
                goto out;
            if (check_sam_write1(conf->out, conf->header, rec, conf->fn_out,
                                 &write_error) < 0)
                goto out;
        }
    }

    if (r < -1) {
        print_error_errno("view", "error reading file \"%s\"", conf->fn_in);
        goto out;
    }

    retval = EXIT_SUCCESS;

 out:
    hts_itr_multi_destroy(iter);
    hts_idx_destroy(conf->hts_idx); // destroy the BAM index
    conf->hts_idx = NULL;
    if (mate_names) {
        // free khash keys
        for (k = 0; k < kh_end(mate_names); ++k)
            if ( kh_exist(mate_names,k) ) free((char*)kh_key(mate_names, k));
        kh_destroy(names,mate_names);
    }
    bam_destroy1(rec);
    return retval;
}

// Common code for processing and writing a record
static inline int process_one_record(samview_settings_t *conf, bam1_t *b,
                                     int *write_error) {
    if (conf->sanitize)
        if (bam_sanitize(conf->header, b, conf->sanitize) < 0)
            return -1;

    if (!process_aln(conf->header, b, conf)) {
        if (!conf->is_count) {
            change_flag(b, conf);
            if (adjust_tags(conf->header, b, conf) != 0)
                return -1;
            if (check_sam_write1(conf->out, conf->header,
                                 b, conf->fn_out, write_error) < 0) {
                return -1;
            }
        }
        conf->count++;
    } else if (conf->unmap) {
        b->core.flag |= BAM_FUNMAP;
        b->core.qual = 0;
        b->core.isize = 0;

        // remove CIGAR
        if (b->core.n_cigar) {
            memmove(bam_get_cigar(b), bam_get_seq(b),
                    b->data + b->l_data - bam_get_seq(b));
            b->l_data -= 4*b->core.n_cigar;
            b->core.n_cigar = 0;
        }

        if (check_sam_write1(conf->out, conf->header,
                             b, conf->fn_out, write_error) < 0) {
            return -1;
        }
    } else {
        if (conf->un_out) {
            if (check_sam_write1(conf->un_out, conf->header,
                                 b, conf->fn_un_out, write_error) < 0) {
                return -1;
            }
        }
    }
    return 0;
}

static int stream_view(samview_settings_t *conf) {
    bam1_t *b = bam_init1();
    int write_error = 0, r;
    if (!b) {
        print_error_errno("view", "could not allocate bam record");
        return 1;
    }
    errno = 0; // prevent false error messages.
    while ((r = sam_read1(conf->in, conf->header, b)) >= 0) {
        if (process_one_record(conf, b, &write_error) < 0) break;
    }
    bam_destroy1(b);
    if (r < -1) {
        print_error_errno("view", "error reading file \"%s\"", conf->fn_in);
        return 1;
    }
    return write_error;
}

static int multi_region_view(samview_settings_t *conf, hts_itr_multi_t *iter)
{
    bam1_t *b = bam_init1();
    int write_error = 0, result;
    if (!b) {
        print_error_errno("view", "could not allocate bam record");
        return 1;
    }
    // fetch alignments
    while ((result = sam_itr_multi_next(conf->in, iter, b)) >= 0) {
        if (process_one_record(conf, b, &write_error) < 0) break;
    }
    hts_itr_multi_destroy(iter);
    bam_destroy1(b);

    if (result < -1) {
        print_error("view", "retrieval of region #%d failed", iter->curr_tid);
        return 1;
    }
    return write_error;
}

// Make mnemonic distinct values for longoption-only options
#define LONGOPT(c)  ((c) + 128)

// Check for ".sam" filenames as sam_open_mode cannot distinguish between
// foo.sam and foo.unknown, both getting mode "".
static int is_sam(const char *fn) {
    if (!fn)
        return 0;
    size_t l = strlen(fn);
    return (l >= 4 && strcasecmp(fn + l-4, ".sam") == 0);
}

static void aux_list_free(samview_settings_t *settings) {
    if (settings->keep_tag)
        kh_destroy(aux_exists, settings->keep_tag);
    if (settings->remove_tag)
        kh_destroy(aux_exists, settings->remove_tag);
}

int main_samview(int argc, char *argv[])
{
    samview_settings_t settings;
    int c, is_header = 0, is_header_only = 0, ret = 0, compress_level = -1, has_index_file = 0, no_pg = 0;
    FILE *fp_out = NULL;
    char out_mode[6] = {0}, out_un_mode[6] = {0};
    char *out_format = "";
    char *arg_list = NULL;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    htsThreadPool p = {NULL, 0};

    memset(&settings,0,sizeof(settings));
    settings.subsam_frac = -1.0;
    settings.count_rf = SAM_FLAG; // don't want 0, and this is quick

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 'T', '@'),
        {"add-flags", required_argument, NULL, LONGOPT('a')},
        {"bam", no_argument, NULL, 'b'},
        {"count", no_argument, NULL, 'c'},
        {"cram", no_argument, NULL, 'C'},
        {"customised-index", no_argument, NULL, 'X'},
        {"customized-index", no_argument, NULL, 'X'},
        {"excl-flags", required_argument, NULL, 'F'},
        {"exclude-flags", required_argument, NULL, 'F'},
        {"expr", required_argument, NULL, 'e'},
        {"expression", required_argument, NULL, 'e'},
        {"fai-reference", required_argument, NULL, 't'},
        {"fast", no_argument, NULL, '1'},
        {"fetch-pairs", no_argument, NULL, 'P'},
        {"header-only", no_argument, NULL, 'H'},
        {"help", no_argument, NULL, LONGOPT('?')},
        {"incl-flags", required_argument, NULL, LONGOPT('g')},
        {"include-flags", required_argument, NULL, LONGOPT('g')},
        {"rf", required_argument, NULL, LONGOPT('g')}, // aka incl-flags
        {"keep-tag", required_argument, NULL, LONGOPT('x') },
        {"library", required_argument, NULL, 'l'},
        {"min-mapq", required_argument, NULL, 'q'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-mq", required_argument, NULL, 'q'},
        {"min-qlen", required_argument, NULL, 'm'},
        {"no-header", no_argument, NULL, LONGOPT('H')},
        {"no-PG", no_argument, NULL, LONGOPT('P')},
        {"output", required_argument, NULL, 'o'},
        {"output-unselected", required_argument, NULL, 'U'},
        {"QNAME-file", required_argument, NULL, 'N'},
        {"qname-file", required_argument, NULL, 'N'},
        {"read-group", required_argument, NULL, 'r'},
        {"read-group-file", required_argument, NULL, 'R'},
        {"readgroup", required_argument, NULL, 'r'},
        {"readgroup-file", required_argument, NULL, 'R'},
        {"region-file", required_argument, NULL, LONGOPT('L')},
        {"regions-file", required_argument, NULL, LONGOPT('L')},
        {"remove-B", no_argument, NULL, 'B'},
        {"remove-flags", required_argument, NULL, LONGOPT('r')},
        {"remove-tag", required_argument, NULL, 'x'},
        {"require-flags", required_argument, NULL, 'f'},
        {"subsample", required_argument, NULL, LONGOPT('s')},
        {"subsample-seed", required_argument, NULL, LONGOPT('S')},
        {"tag", required_argument, NULL, 'd'},
        {"tag-file", required_argument, NULL, 'D'},
        {"target-file", required_argument, NULL, 'L'},
        {"targets-file", required_argument, NULL, 'L'},
        {"uncompressed", no_argument, NULL, 'u'},
        {"unmap", no_argument, NULL, 'p'},
        {"unoutput", required_argument, NULL, 'U'},
        {"use-index", no_argument, NULL, 'M'},
        {"with-header", no_argument, NULL, 'h'},
        {"sanitize", required_argument, NULL, 'z'},
    };

    /* parse command-line options */
    strcpy(out_mode, "w");
    strcpy(out_un_mode, "w");
    if (argc == 1 && isatty(STDIN_FILENO))
        return usage(stdout, EXIT_SUCCESS, 0);

    // Suppress complaints about '?' being an unrecognised option.  Without
    // this we have to put '?' in the options list, which makes it hard to
    // tell a bad long option from the use of '-?' (both return '?' and
    // set optopt to '\0').
    opterr = 0;

    char *tmp;
    while ((c = getopt_long(argc, argv,
                            "SbBcCt:h1Ho:O:q:f:F:G:ul:r:T:R:N:d:D:L:s:@:m:x:U:MXe:pPz:",
                            lopts, NULL)) >= 0) {
        switch (c) {
        case 's':
            settings.subsam_seed = strtol(optarg, &tmp, 10);
            if (tmp && *tmp == '.') {
                settings.subsam_frac = strtod(tmp, &tmp);
                if (*tmp) ret = 1;
            } else {
                ret = 1;
            }

            if (ret == 1) {
                print_error("view", "Incorrect sampling argument \"%s\"", optarg);
                goto view_end;
            }
            settings.count_rf |= SAM_QNAME;
            break;
        case LONGOPT('s'):
            settings.subsam_frac = strtod(optarg, &tmp);
            if (*tmp || settings.subsam_frac < 0.0 || settings.subsam_frac > 1.0) {
                print_error("view", "Incorrect sampling argument \"%s\"", optarg);
                goto view_end;
            }
            settings.count_rf |= SAM_QNAME;
            break;
        case LONGOPT('S'): settings.subsam_seed = atoi(optarg); break;
        case 'm':
            settings.min_qlen = atoi(optarg);
            settings.count_rf |= SAM_SEQ;
            break;
        case 'c': settings.is_count = 1; break;
        case 'S': break;
        case 'b': out_format = "b"; break;
        case 'C': out_format = "c"; break;
        case 't': settings.fn_fai = strdup(optarg); break;
        case 'h': is_header = 1; break;
        case 'H': is_header_only = 1; break;
        case LONGOPT('H'): is_header = is_header_only = 0; break;
        case 'o': settings.fn_out = strdup(optarg); break;
        case 'U': settings.fn_un_out = strdup(optarg); break;
        case 'X': has_index_file = 1; break;
        case 'f':
            settings.flag_on |= bam_str2flag(optarg);
            settings.count_rf |= SAM_FLAG | SAM_RNEXT;
            break;
        case 'F':
            settings.flag_off |= bam_str2flag(optarg);
            settings.count_rf |= SAM_FLAG | SAM_RNEXT;
            break;
        case LONGOPT('g'):
            settings.flag_anyon |= bam_str2flag(optarg);
            settings.count_rf |= SAM_FLAG | SAM_RNEXT;
            break;
        case 'G':
            settings.flag_alloff |= bam_str2flag(optarg);
            settings.count_rf |= SAM_FLAG | SAM_RNEXT;
            break;
        case 'q':
            settings.min_mapQ = atoi(optarg);
            settings.count_rf |= SAM_MAPQ;
            break;
        case 'u': compress_level = 0; break;
        case '1': compress_level = 1; break;
        case 'l':
            settings.library = strdup(optarg);
            settings.count_rf |= SAM_RGAUX;
            break;
        case 'p': settings.unmap = 1; break;
        case 'P': settings.fetch_pairs = 1; settings.multi_region = 1; break;
        case 'z':
            if ((settings.sanitize = bam_sanitize_options(optarg)) < 0) {
                ret = 1;
                goto view_end;
            }
            break;
        case LONGOPT('L'):
            settings.multi_region = 1;
            // fall through
        case 'L':
            if ((settings.bed = bed_read(optarg)) == NULL) {
                print_error_errno("view", "Could not read file \"%s\"", optarg);
                ret = 1;
                goto view_end;
            }
            settings.count_rf |= SAM_POS | SAM_RNAME | SAM_CIGAR;
            break;
        case 'r':
            if (add_read_group_single("view", &settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            settings.count_rf |= SAM_RGAUX;
            break;
        case 'R':
            if (add_read_groups_file("view", &settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            settings.count_rf |= SAM_RGAUX;
            break;
        case 'N':
            if (add_read_names_file("view", &settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            settings.count_rf |= SAM_QNAME;
            break;

        case 'd':
            if (strlen(optarg) < 2 || (strlen(optarg) > 2 && optarg[2] != ':')) {
                print_error("view", "Invalid \"tag:value\" option: \"%s\"", optarg);
                ret = 1;
                goto view_end;
            }

            if (settings.tag) {
                if (settings.tag[0] != optarg[0] || settings.tag[1] != optarg[1]) {
                    print_error("view", "Different tag \"%s\" was specified before: \"%s\"", settings.tag, optarg);
                    ret = 1;
                    goto view_end;
                }
            } else {
                if (!(settings.tag = calloc(3, 1))) {
                    print_error("view", "Could not allocate memory for tag: \"%s\"", optarg);
                    ret = 1;
                    goto view_end;
                }
                memcpy(settings.tag, optarg, 2);
            }

            if (strlen(optarg) > 3 && add_tag_value_single("view", &settings, optarg+3) != 0) {
                print_error("view", "Could not add tag:value \"%s\"", optarg);
                ret = 1;
                goto view_end;
            }
            // Some tag filtering affects other fields
            if (memcmp(settings.tag, "NM", 2) == 0 ||
                memcmp(settings.tag, "MD", 2) == 0)
                settings.count_rf |= SAM_AUX | SAM_SEQ;
            else if (memcmp(settings.tag, "RG", 2) == 0)
                settings.count_rf |= SAM_RGAUX;
            else
                settings.count_rf |= SAM_AUX;
            break;

        case 'D':
            // Allow ";" as delimiter besides ":" to support MinGW CLI POSIX
            // path translation as described at:
            // http://www.mingw.org/wiki/Posix_path_conversion
            if (strlen(optarg) < 4 || (optarg[2] != ':' && optarg[2] != ';')) {
                print_error("view", "Invalid \"tag:file\" option: \"%s\"", optarg);
                ret = 1;
                goto view_end;
            }

            if (settings.tag) {
                if (settings.tag[0] != optarg[0] || settings.tag[1] != optarg[1]) {
                    print_error("view", "Different tag \"%s\" was specified before: \"%s\"", settings.tag, optarg);
                    ret = 1;
                    goto view_end;
                }
            } else {
                if (!(settings.tag = calloc(3, 1))) {
                    print_error("view", "Could not allocate memory for tag: \"%s\"", optarg);
                    ret = 1;
                    goto view_end;
                }
                memcpy(settings.tag, optarg, 2);
            }

            if (add_tag_values_file("view", &settings, optarg+3) != 0) {
                ret = 1;
                goto view_end;
            }
            // Some tag filtering affects other fields
            if (memcmp(settings.tag, "NM", 2) == 0 ||
                memcmp(settings.tag, "MD", 2) == 0)
                settings.count_rf |= SAM_AUX | SAM_SEQ;
            else if (memcmp(settings.tag, "RG", 2) == 0)
                settings.count_rf |= SAM_RGAUX;
            else
                settings.count_rf |= SAM_AUX;
            break;

        case LONGOPT('?'):
            return usage(stdout, EXIT_SUCCESS, 1);
        case '?':
            if (optopt == '?') {  // '-?' appeared on command line
                return usage(stdout, EXIT_SUCCESS, 1);
            } else {
                if (optopt) { // Bad short option
                    print_error("view", "invalid option -- '%c'", optopt);
                } else { // Bad long option
                    // Do our best.  There is no good solution to finding
                    // out what the bad option was.
                    // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
                    if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
                        print_error("view", "unrecognised option '%s'",
                                    argv[optind - 1]);
                    }
                }
                return usage(stderr, EXIT_FAILURE, 0);
            }
        case 'B': settings.remove_B = 1; break;

        case 'M': settings.multi_region = 1; break;
        case LONGOPT('P'): no_pg = 1; break;
        case 'e':
            if (!(settings.filter = hts_filter_init(optarg))) {
                print_error("main_samview", "Couldn't initialise filter");
                return 1;
            }
            settings.count_rf = INT_MAX; // no way to know what we need
            break;
        case LONGOPT('r'): settings.remove_flag |= bam_str2flag(optarg); break;
        case LONGOPT('a'): settings.add_flag |= bam_str2flag(optarg); break;

        case 'x':
            if (*optarg == '^') {
                if (parse_aux_list(&settings.keep_tag, optarg+1, "main_samview")) {
                    aux_list_free(&settings);
                    return usage(stderr, EXIT_FAILURE, 0);
                }
            } else {
                if (parse_aux_list(&settings.remove_tag, optarg, "main_samview")) {
                    aux_list_free(&settings);
                    return usage(stderr, EXIT_FAILURE, 0);
                }
            }
            break;

        case LONGOPT('x'):
            if (parse_aux_list(&settings.keep_tag, optarg, "main_samview")) {
                aux_list_free(&settings);
                return usage(stderr, EXIT_FAILURE, 0);
            }
            break;

        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) != 0)
                return usage(stderr, EXIT_FAILURE, 0);
            break;
        }
    }
    if (settings.is_count && settings.fetch_pairs)
    {
        print_error("view","The options -P and -c cannot be combined\n");
        return 1;
    }
    if (settings.fn_fai == 0 && ga.reference) settings.fn_fai = fai_path(ga.reference);
    if (is_header_only) is_header = 1;
    // File format auto-detection first
    if (settings.fn_out)    sam_open_mode(out_mode+1,    settings.fn_out,    NULL);
    if (settings.fn_un_out) sam_open_mode(out_un_mode+1, settings.fn_un_out, NULL);

    // -1 or -u without an explicit format (-b, -C) => check fn extensions
    if (!*out_format && compress_level >= 0) {
        if (compress_level == 0 &&
            (out_mode[strlen(out_mode)-1] == 'z' ||
             out_un_mode[strlen(out_un_mode)-1] == 'z'))
            // z, fz, Fz sanity check
            fprintf(stderr, "[view] Warning option -u ignored due to"
                    " filename suffix\n");

        // If known extension, use it, otherwise BAM
        if (!(out_mode[1] || is_sam(settings.fn_out)))
            out_mode[1] = 'b';

        if (!(out_un_mode[1] || is_sam(settings.fn_un_out)))
            out_un_mode[1] = 'b';
    } else if (*out_format) {
        out_mode[1] = out_un_mode[1] = *out_format;
    }

    if (compress_level >= 0) {
        char tmp[2];
        tmp[0] = compress_level + '0'; tmp[1] = '\0';
        strcat(out_mode, tmp);
        strcat(out_un_mode, tmp);
    }
    if (argc == optind && isatty(STDIN_FILENO)) {
        print_error("view", "No input provided or missing option argument.");
        return usage(stderr, EXIT_FAILURE, 0); // potential memory leak...
    }

    if (settings.unmap && settings.fn_un_out) {
        print_error("view", "Options --unoutput and --unmap are mutually exclusive.");
        ret = 1;
        goto view_end;
    }

    if (settings.subsam_seed != 0) {
        // Convert likely user input 1,2,... to pseudo-random
        // values with more entropy and more bits set
        srand(settings.subsam_seed);
        settings.subsam_seed = rand();
    }

    settings.fn_in = (optind < argc)? argv[optind] : "-";
    if ((settings.in = sam_open_format(settings.fn_in, "r", &ga.in)) == 0) {
        print_error_errno("view", "failed to open \"%s\" for reading", settings.fn_in);
        ret = 1;
        goto view_end;
    }

    if (settings.fn_fai) {
        if (hts_set_fai_filename(settings.in, settings.fn_fai) != 0) {
            fprintf(stderr, "[main_samview] failed to use reference \"%s\".\n", settings.fn_fai);
            ret = 1;
            goto view_end;
        }
    }
    if ((settings.header = sam_hdr_read(settings.in)) == 0) {
        fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", settings.fn_in);
        ret = 1;
        goto view_end;
    }
    if (settings.rghash) {
        sam_hdr_remove_lines(settings.header, "RG", "ID", settings.rghash);
    }
    if (!settings.is_count) {
        if ((settings.out = sam_open_format(settings.fn_out? settings.fn_out : "-", out_mode, &ga.out)) == 0) {
            print_error_errno("view", "failed to open \"%s\" for writing", settings.fn_out? settings.fn_out : "standard output");
            ret = 1;
            goto view_end;
        }
        if (settings.fn_fai) {
            if (hts_set_fai_filename(settings.out, settings.fn_fai) != 0) {
                fprintf(stderr, "[main_samview] failed to use reference \"%s\".\n", settings.fn_fai);
                ret = 1;
                goto view_end;
            }
        }
        autoflush_if_stdout(settings.out, settings.fn_out);

        if (!no_pg) {
            if (!(arg_list = stringify_argv(argc+1, argv-1))) {
                print_error("view", "failed to create arg_list");
                ret = 1;
                goto view_end;
            }
            if (sam_hdr_add_pg(settings.header, "samtools",
                                         "VN", samtools_version(),
                                         arg_list ? "CL": NULL,
                                         arg_list ? arg_list : NULL,
                                         NULL)) {
                print_error("view", "failed to add PG line to the header");
                ret = 1;
                goto view_end;
            }
        }

        if (ga.write_index || is_header ||
            out_mode[1] == 'b' || out_mode[1] == 'c' ||
            (ga.out.format != sam && ga.out.format != unknown_format))  {
            if (sam_hdr_write(settings.out, settings.header) != 0) {
                fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                ret = 1;
                goto view_end;
            }
        }
        if (ga.write_index) {
            if (!(settings.fn_out_idx = auto_index(settings.out, settings.fn_out, settings.header))) {
                ret = 1;
                goto view_end;
            }
        }

        if (settings.fn_un_out) {
            if ((settings.un_out = sam_open_format(settings.fn_un_out, out_un_mode, &ga.out)) == 0) {
                print_error_errno("view", "failed to open \"%s\" for writing", settings.fn_un_out);
                ret = 1;
                goto view_end;
            }
            if (settings.fn_fai) {
                if (hts_set_fai_filename(settings.un_out, settings.fn_fai) != 0) {
                    fprintf(stderr, "[main_samview] failed to use reference \"%s\".\n", settings.fn_fai);
                    ret = 1;
                    goto view_end;
                }
            }
            autoflush_if_stdout(settings.un_out, settings.fn_un_out);
            if (ga.write_index || is_header ||
                out_un_mode[1] == 'b' || out_un_mode[1] == 'c' ||
                (ga.out.format != sam && ga.out.format != unknown_format))  {
                if (sam_hdr_write(settings.un_out, settings.header) != 0) {
                    fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                    ret = 1;
                    goto view_end;
                }
            }
            if (ga.write_index) {
                if (!(settings.fn_un_out_idx = auto_index(settings.un_out, settings.fn_un_out, settings.header))) {
                    ret = 1;
                    goto view_end;
                }
            }
        }
    }
    else {
        if (settings.fn_out) {
            fp_out = fopen(settings.fn_out, "w");
            if (fp_out == NULL) {
                print_error_errno("view", "can't create \"%s\"", settings.fn_out);
                ret = EXIT_FAILURE;
                goto view_end;
            }
        }
        settings.unmap = 0;  // Not valid in counting mode
    }

    if (ga.nthreads > 0) {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            ret = 1;
            goto view_end;
        }
        hts_set_opt(settings.in,  HTS_OPT_THREAD_POOL, &p);
        if (settings.out) hts_set_opt(settings.out, HTS_OPT_THREAD_POOL, &p);
    }
    if (is_header_only) goto view_end; // no need to print alignments


    // Initialize BAM/CRAM index
    char **regs = NULL;
    int nregs = 0;
    if ( has_index_file && optind <= argc - 2 ) {
        regs = optind < argc-2 ? &argv[optind+2] : NULL;
        nregs = argc - optind - 2;
        settings.fn_idx_in = argv[optind+1];
    } else if (!has_index_file && optind < argc - 1 ) {
        regs = &argv[optind+1];
        nregs = argc - optind - 1;
    } else if ( has_index_file && argc-optind < 2) {
        print_error("view", "Incorrect number of arguments for -X option. Aborting.");
        return 1;
    }
    if (regs)
        settings.count_rf |= SAM_POS | SAM_RNAME | SAM_CIGAR;

    if ( settings.fn_idx_in || nregs || settings.multi_region )
    {
        settings.hts_idx = settings.fn_idx_in ? sam_index_load2(settings.in, settings.fn_in, settings.fn_idx_in) : sam_index_load(settings.in, settings.fn_in);
        if ( !settings.hts_idx )
        {
            print_error("view", "Random alignment retrieval only works for indexed SAM.gz, BAM or CRAM files.");
            return 1;
        }
    }

    if (settings.is_count)
        // Won't fail, but also wouldn't matter if it did
        hts_set_opt(settings.in, CRAM_OPT_REQUIRED_FIELDS, settings.count_rf);

    if ( settings.fetch_pairs )
    {
        hts_itr_multi_t *iter = multi_region_init(&settings, regs, nregs);
        ret = iter ? fetch_pairs_collect_mates(&settings, iter) : 1;
        if (ret) goto view_end;
    }
    else if ( settings.multi_region )
    {
        hts_itr_multi_t *iter = multi_region_init(&settings, regs, nregs);
        ret = iter ? multi_region_view(&settings, iter) : 1;
        if (ret) goto view_end;
    }
    else if ( !settings.hts_idx || optind+1 >= argc-has_index_file ) {
        // stream through the entire file
        ret = stream_view(&settings);
        if (ret) goto view_end;
    } else {   // retrieve alignments in specified regions
        int i;
        for (i = (has_index_file)? optind+2 : optind+1; i < argc; ++i) {
            hts_itr_t *iter = sam_itr_querys(settings.hts_idx, settings.header, argv[i]); // parse a region in the format like `chr2:100-200'
            if (iter == NULL) { // region invalid or reference name not found
                fprintf(stderr, "[main_samview] region \"%s\" specifies an invalid region or unknown reference. Continue anyway.\n", argv[i]);
                continue;
            }
            // fetch alignments
            ret = multi_region_view(&settings, iter);
            if (ret) goto view_end;
        }
    }

    if ( settings.hts_idx ) hts_idx_destroy(settings.hts_idx);

    if (ga.write_index) {
        if (sam_idx_save(settings.out) < 0) {
            print_error_errno("view", "writing index failed");
            ret = 1;
        }
        if (settings.un_out && sam_idx_save(settings.un_out) < 0) {
            print_error_errno("view", "writing index failed");
            ret = 1;
        }
    }

view_end:
    if (settings.is_count && ret == 0) {
        if (fprintf(settings.fn_out? fp_out : stdout, "%" PRId64 "\n", settings.count) < 0) {
            if (settings.fn_out) print_error_errno("view", "writing to \"%s\" failed", settings.fn_out);
            else print_error_errno("view", "writing to standard output failed");
            ret = EXIT_FAILURE;
        }
    }

    // close files, free and return
    if (settings.in) check_sam_close("view", settings.in, settings.fn_in, "standard input", &ret);
    if (settings.out) check_sam_close("view", settings.out, settings.fn_out, "standard output", &ret);
    if (settings.un_out) check_sam_close("view", settings.un_out, settings.fn_un_out, "file", &ret);
    if (fp_out) fclose(fp_out);

    free(settings.fn_fai); free(settings.fn_out); free(settings.library);  free(settings.fn_un_out);
    sam_global_args_free(&ga);
    if ( settings.header ) sam_hdr_destroy(settings.header);
    if (settings.bed) bed_destroy(settings.bed);
    if (settings.rghash) {
        khint_t k;
        for (k = 0; k < kh_end(settings.rghash); ++k)
            if (kh_exist(settings.rghash, k)) free((char*)kh_key(settings.rghash, k));
        kh_destroy(str, settings.rghash);
    }
    if (settings.rnhash) {
        khint_t k;
        for (k = 0; k < kh_end(settings.rnhash); ++k)
            if (kh_exist(settings.rnhash, k)) free((char*)kh_key(settings.rnhash, k));
        kh_destroy(str, settings.rnhash);
    }
    if (settings.tvhash) {
        khint_t k;
        for (k = 0; k < kh_end(settings.tvhash); ++k)
            if (kh_exist(settings.tvhash, k)) free((char*)kh_key(settings.tvhash, k));
        kh_destroy(str, settings.tvhash);
    }

    if (settings.remove_aux_len) {
        free(settings.remove_aux);
    }
    if (settings.tag) {
        free(settings.tag);
    }
    if (settings.filter)
        hts_filter_free(settings.filter);

    if (p.pool)
        hts_tpool_destroy(p.pool);

    if (settings.fn_out_idx)
        free(settings.fn_out_idx);
    if (settings.fn_un_out_idx)
        free(settings.fn_un_out_idx);
    free(arg_list);

    aux_list_free(&settings);

    return ret;
}

static int usage(FILE *fp, int exit_status, int is_long_help)
{
    fprintf(fp,
"\n"
"Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]\n"
"\n"

"Output options:\n"
"  -b, --bam                  Output BAM\n"
"  -C, --cram                 Output CRAM (requires -T)\n"
"  -1, --fast                 Use fast BAM compression (and default to --bam)\n"
"  -u, --uncompressed         Uncompressed BAM output (and default to --bam)\n"
"  -h, --with-header          Include header in SAM output\n"
"  -H, --header-only          Print SAM header only (no alignments)\n"
"      --no-header            Print SAM alignment records only [default]\n"
"  -c, --count                Print only the count of matching records\n"
"  -o, --output FILE          Write output to FILE [standard output]\n"
"  -U, --unoutput FILE, --output-unselected FILE\n"
"                             Output reads not selected by filters to FILE\n"
"  -p, --unmap                Set flag to UNMAP on reads not selected\n"
"                             then write to output file.\n"
"  -P, --fetch-pairs          Retrieve complete pairs even when outside of region\n"
"Input options:\n"
"  -t, --fai-reference FILE   FILE listing reference names and lengths\n"
"  -M, --use-index            Use index and multi-region iterator for regions\n"
"      --region[s]-file FILE  Use index to include only reads overlapping FILE\n"
"  -X, --customized-index     Expect extra index file argument after <in.bam>\n"
"\n"
"Filtering options (Only include in output reads that...):\n"
"  -L, --target[s]-file FILE  ...overlap (BED) regions in FILE\n"
"  -r, --read-group STR       ...are in read group STR\n"
"  -R, --read-group-file FILE ...are in a read group listed in FILE\n"
"  -N, --qname-file FILE      ...whose read name is listed in FILE\n"
"  -d, --tag STR1[:STR2]      ...have a tag STR1 (with associated value STR2)\n"
"  -D, --tag-file STR:FILE    ...have a tag STR whose value is listed in FILE\n"
"  -q, --min-MQ INT           ...have mapping quality >= INT\n"
"  -l, --library STR          ...are in library STR\n"
"  -m, --min-qlen INT         ...cover >= INT query bases (as measured via CIGAR)\n"
"  -e, --expr STR             ...match the filter expression STR\n"
"  -f, --require-flags FLAG   ...have all of the FLAGs present\n"             //   F&x == x
"  -F, --excl[ude]-flags FLAG ...have none of the FLAGs present\n"            //   F&x == 0
"      --rf, --incl-flags, --include-flags FLAG\n"
"                             ...have some of the FLAGs present\n"
"  -G FLAG                    EXCLUDE reads with all of the FLAGs present\n"  // !(F&x == x)  TODO long option
"      --subsample FLOAT      Keep only FLOAT fraction of templates/read pairs\n"
"      --subsample-seed INT   Influence WHICH reads are kept in subsampling [0]\n"
"  -s INT.FRAC                Same as --subsample 0.FRAC --subsample-seed INT\n"
"\n"
"Processing options:\n"
"      --add-flags FLAG       Add FLAGs to reads\n"
"      --remove-flags FLAG    Remove FLAGs from reads\n"
"  -x, --remove-tag STR\n"
"               Comma-separated read tags to strip (repeatable) [null]\n"
"      --keep-tag STR\n"
"               Comma-separated read tags to preserve (repeatable) [null].\n"
"               Equivalent to \"-x ^STR\"\n"
"  -B, --remove-B             Collapse the backward CIGAR operation\n"
"  -z, --sanitize FLAGS       Perform sanitity checking and fixing on records.\n"
"                             FLAGS is comma separated (see manual). [off]\n"
"\n"
"General options:\n"
"  -?, --help   Print long help, including note about region specification\n"
"  -S           Ignored (input format is auto-detected)\n"
"      --no-PG  Do not add a PG line\n");

    sam_global_opt_help(fp, "-.O.T@..");
    fprintf(fp, "\n");

    if (is_long_help)
        fprintf(fp,
"Notes:\n"
"\n"
"1. This command now auto-detects the input format (BAM/CRAM/SAM).\n"
"   Further control over the CRAM format can be specified by using the\n"
"   --output-fmt-option, e.g. to specify the number of sequences per slice\n"
"   and to use avoid reference based compression:\n"
"\n"
"\tsamtools view -C --output-fmt-option seqs_per_slice=5000 \\\n"
"\t   --output-fmt-option no_ref -o out.cram in.bam\n"
"\n"
"   Options can also be specified as a comma separated list within the\n"
"   --output-fmt value too.  For example this is equivalent to the above\n"
"\n"
"\tsamtools view --output-fmt cram,seqs_per_slice=5000,no_ref \\\n"
"\t   -o out.cram in.bam\n"
"\n"
"2. The file supplied with `-t' is SPACE/TAB delimited with the first\n"
"   two fields of each line consisting of the reference name and the\n"
"   corresponding sequence length. The `.fai' file generated by \n"
"   `samtools faidx' is suitable for use as this file. This may be an\n"
"   empty file if reads are unaligned.\n"
"\n"
"3. SAM->BAM conversion:  samtools view -bT ref.fa in.sam.gz\n"
"\n"
"4. BAM->SAM conversion:  samtools view -h in.bam\n"
"\n"
"5. A region should be presented in one of the following formats:\n"
"   `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n"
"   specified, the input alignment file must be a sorted and indexed\n"
"   alignment (BAM/CRAM) file.\n"
"\n"
"6. Option `-u' is preferred over `-b' when the output is piped to\n"
"   another samtools command.\n"
"\n"
"7. Option `-M`/`--use-index` causes overlaps with `-L` BED file regions and\n"
"   command-line region arguments to be computed using the multi-region iterator\n"
"   and an index. This increases speed, omits duplicates, and outputs the reads\n"
"   as they are ordered in the input SAM/BAM/CRAM file.\n"
"\n"
"8. Options `-L`/`--target[s]-file` and `--region[s]-file` may not be used\n"
"   together. `--region[s]-file FILE` is simply equivalent to `-M -L FILE`,\n"
"   so using both causes one of the specified BED files to be ignored.\n"
"\n");

    return exit_status;
}

static int head_usage(FILE *fp, int exit_status)
{
    fprintf(fp,
"Usage: samtools head [OPTION]... [FILE]\n"
"Options:\n"
"  -h, --headers INT   Display INT header lines [all]\n"
"  -n, --records INT   Display INT alignment record lines [none]\n"
);
    sam_global_opt_help(fp, "-.--T@-.");
    return exit_status;
}

int main_head(int argc, char *argv[])
{
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 'T', '@'),
        { "headers", required_argument, NULL, 'h' },
        { "records", required_argument, NULL, 'n' },
        { NULL, 0, NULL, 0 }
    };
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    int all_headers = 1;
    uint64_t nheaders = 0;
    uint64_t nrecords = 0;

    int c, nargs;
    while ((c = getopt_long(argc, argv, "h:n:T:@:", lopts, NULL)) >= 0)
        switch (c) {
        case 'h': all_headers = 0; nheaders = strtoull(optarg, NULL, 0); break;
        case 'n': nrecords = strtoull(optarg, NULL, 0); break;
        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            return head_usage(stderr, EXIT_FAILURE);
        }

    nargs = argc - optind;
    if (nargs == 0 && isatty(STDIN_FILENO))
        return head_usage(stdout, EXIT_SUCCESS);
    else if (nargs > 1)
        return head_usage(stderr, EXIT_FAILURE);

    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    kstring_t str = KS_INITIALIZE;
    bam1_t *b = NULL;

    const char *fname = (nargs == 1)? argv[optind] : "-";
    fp = sam_open_format(fname, "r", &ga.in);
    if (fp == NULL) {
        if (strcmp(fname, "-") != 0)
            print_error_errno("head", "failed to open \"%s\" for reading", fname);
        else
            print_error_errno("head", "failed to open standard input for reading");
        goto err;
    }

    if (ga.nthreads > 0) hts_set_threads(fp, ga.nthreads);

    hdr = sam_hdr_read(fp);
    if (hdr == NULL) {
        if (strcmp(fname, "-") != 0)
            print_error("head", "failed to read the header from \"%s\"", fname);
        else
            print_error("head", "failed to read the header");
        goto err;
    }

    if (all_headers) {
        fputs(sam_hdr_str(hdr), stdout);
    }
    else if (nheaders > 0) {
        const char *text = sam_hdr_str(hdr);
        const char *lim = text;
        uint64_t n;
        for (n = 0; n < nheaders; n++) {
            lim = strchr(lim, '\n');
            if (lim) lim++;
            else break;
        }
        if (lim) fwrite(text, lim - text, 1, stdout);
        else fputs(text, stdout);
    }

    if (nrecords > 0) {
        b = bam_init1();
        uint64_t n;
        int r;
        for (n = 0; n < nrecords && (r = sam_read1(fp, hdr, b)) >= 0; n++) {
            if (sam_format1(hdr, b, &str) < 0) {
                print_error_errno("head", "couldn't format record");
                goto err;
            }
            puts(ks_str(&str));
        }
        if (r < -1) {
            print_error("head", "\"%s\" is truncated", fname);
            goto err;
        }
        bam_destroy1(b);
        ks_free(&str);
    }

    sam_hdr_destroy(hdr);
    sam_close(fp);
    sam_global_args_free(&ga);

    return EXIT_SUCCESS;

err:
    if (fp) sam_close(fp);
    sam_hdr_destroy(hdr);
    bam_destroy1(b);
    ks_free(&str);
    sam_global_args_free(&ga);
    return EXIT_FAILURE;
}
