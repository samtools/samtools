/*  stats.c -- This is the former bamcheck integrated into samtools/htslib.

    Copyright (C) 2020-2021 Genome Research Ltd.

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

/*
 * This tool is designed to give "samtools stats" style output, but dedicated
 * to small amplicon sequencing projects.  It gathers stats on the
 * distribution of reads across amplicons.
 */

/*
 * TODO:
 * - Cope with multiple references.  What do we do here?  Just request one?
 * - Permit regions rather than consuming whole file (maybe solves above).
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <unistd.h>
#include <math.h>

#include <htslib/sam.h>
#include <htslib/khash.h>

#include "samtools.h"
#include "sam_opts.h"
#include "bam_ampliconclip.h"

KHASH_MAP_INIT_INT64(tcoord, int64_t)
KHASH_MAP_INIT_STR(qname, int64_t)

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(a) ((a)>=0?(a):-(a))
#endif

#define TCOORD_MIN_COUNT   10
#define MAX_AMP 1000       // Default maximum number of amplicons
#define MAX_AMP_LEN 1000   // Default maximum length of any single amplicon
#define MAX_PRIMER_PER_AMPLICON 4  // Max primers per LEFT/RIGHT
#define MAX_DEPTH 5        // Number of different depths permitted

typedef struct {
    sam_global_args ga;
    uint32_t flag_require;
    uint32_t flag_filter;
    int max_delta;   // Used for matching read to amplicon primer loc
    int min_depth[MAX_DEPTH]; // Used for coverage; must be >= min_depth deep
    int use_sample_name;
    int max_amp;     // Total number of amplicons
    int max_amp_len; // Maximum length of an individual amplicon
    double depth_bin;// aggregate depth within this fraction
    int tlen_adj;    // Adjust tlen by this amount, due to clip but no fixmate
    FILE *out_fp;
    char *argv;
    int tcoord_min_count;
    int tcoord_bin;
    int multi_ref;
} astats_args_t;

typedef struct {
    int nseq;       // total sequence count
    int nfiltered;  // sequence filtered
    int nfailprimer;// count of sequences not matching the primer locations

    // Sizes of memory allocated below, to permit reset
    int max_amp, max_amp_len, max_len;

    // Summary across all samples, sum(x) plus sum(x^2) for s.d. calc
    int64_t *nreads, *nreads2;          // [max_amp]
    double  *nfull_reads;               // [max_amp]; 0.5/read if paired.
    double  *nrperc, *nrperc2;          // [max_amp]
    int64_t *nbases, *nbases2;          // [max_amp]
    int64_t *coverage;                  // [max_amp][max_amp_len]
    double  (*covered_perc)[MAX_DEPTH]; // [max_amp][MAX_DEPTH]
    double  (*covered_perc2)[MAX_DEPTH];// [max_amp][MAX_DEPTH];
    khash_t(tcoord) **tcoord;           // [max_amp+1]

    // 0 is correct pair, 1 is incorrect pair, 2 is unidentified
    int     (*amp_dist)[3];             // [MAX_AMP][3];

    int *depth_valid; // [max_len]
    int *depth_all;   // [max_len]
    khash_t(qname) *qend;  // queryname end, for overlap removal
} astats_t;

// We can have multiple primers for LEFT / RIGHT, so this
// permits detection by any compatible combination.
// One reference:
typedef struct {
    int64_t left[MAX_PRIMER_PER_AMPLICON];
    int nleft;
    int64_t right[MAX_PRIMER_PER_AMPLICON];
    int nright;
    int64_t max_left, min_right; // inner dimensions
    int64_t min_left, max_right; // outer dimensions
} amplicon_t;

// Multiple references, we have an array of amplicons_t - one per used ref.
// We have per reference local and global stats here, as some of the stats
// are coordinate based.  However we report them combined together as a single
// list across all references.
// "namp" is the number of amplicons in this reference, but they're
// numbered first_amp to first_amp+namp-1 inclusively.
typedef struct {
    int tid, namp;
    int64_t len;
    bed_entry_list_t *sites;
    amplicon_t *amp;
    astats_t *lstats, *gstats; // local (1 file) and global (all file) stats
    const char *ref;           // ref name (pointer to the bed hash table key)
    int first_amp;             // first amplicon number for this ref
} amplicons_t;

// Reinitialised for each new reference/chromosome.
// Counts from 1 to namp, -1 for no match and 0 for ?.
static int *pos2start = NULL;
static int *pos2end = NULL;
static int pos2size = 0; // allocated size of pos2start/end

// Lookup table to go from position to amplicon based on
// read start / end.
static int initialise_amp_pos_lookup(astats_args_t *args,
                                     amplicons_t *amps,
                                     int ref) {
    int64_t i, j;
    amplicon_t *amp = amps[ref].amp;
    int64_t max_len = amps[ref].len;
    int namp = amps[ref].namp;

    if (max_len+1 > pos2size) {
        if (!(pos2start = realloc(pos2start, (max_len+1)*sizeof(*pos2start))))
            return -1;
        if (!(pos2end   = realloc(pos2end,   (max_len+1)*sizeof(*pos2end))))
            return -1;
        pos2size = max_len;
    }
    for (i = 0; i < max_len; i++)
        pos2start[i] = pos2end[i] = -1;

    for (i = 0; i < namp; i++) {
        for (j = 0; j < amp[i].nleft; j++) {
            int64_t p;
            for (p = amp[i].left[j] - args->max_delta;
                 p <= amp[i].left[j] + args->max_delta; p++) {
                if (p < 1 || p > max_len)
                    continue;
                pos2start[p-1] = i;
            }
        }
        for (j = 0; j < amp[i].nright; j++) {
            int64_t p;
            for (p = amp[i].right[j] - args->max_delta;
                 p <= amp[i].right[j] + args->max_delta; p++) {
                if (p < 1 || p > max_len)
                    continue;
                pos2end[p-1] = i;
            }
        }
    }

    return 0;
}

// Counts amplicons.
// Assumption: input BED file alternates between LEFT and RIGHT primers
// per amplicon, thus we can count the number based on the switching
// orientation.
static int count_amplicon(bed_entry_list_t *sites) {
    int i, namp, last_rev = 0;
    for (i = namp = 0; i < sites->length; i++) {
        if (sites->bp[i].rev == 0 && last_rev)
            namp++;
        last_rev = sites->bp[i].rev;
    }

    return ++namp;
}

// We're only interest in the internal part of the amplicon.
// Our bed file has LEFT start/end followed by RIGHT start/end,
// so collapse these to LEFT end / RIGHT start.
//
// Returns right most amplicon position on success,
//         < 0 on error
static int64_t bed2amplicon(astats_args_t *args, bed_entry_list_t *sites,
                            amplicon_t *amp, int *namp, int do_title,
                            const char *ref, int first_amp) {
    int i, j;
    int64_t max_right = 0;
    FILE *ofp = args->out_fp;

    *namp = 0;

    // Assume all primers for the same amplicon are adjacent in BED
    // with all + followed by all -.  Thus - to + signifies next primer set.
    int last_rev = 0;
    amp[0].max_left = 0;
    amp[0].min_right = INT64_MAX;
    amp[0].min_left = INT64_MAX;
    amp[0].max_right = 0;
    if (do_title) {
        fprintf(ofp, "# Amplicon locations from BED file.\n");
        fprintf(ofp, "# LEFT/RIGHT are <start>-<end> format and "
                "comma-separated for alt-primers.\n");
        if (args->multi_ref)
            fprintf(ofp, "#\n# AMPLICON\tREF\tNUMBER\tLEFT\tRIGHT\n");
        else
            fprintf(ofp, "#\n# AMPLICON\tNUMBER\tLEFT\tRIGHT\n");
    }
    for (i = j = 0; i < sites->length; i++) {
        if (i == 0 && sites->bp[i].rev != 0) {
            fprintf(stderr, "[ampliconstats] error: BED file should start"
                    " with the + strand primer\n");
            return -1;
        }
        if (sites->bp[i].rev == 0 && last_rev) {
            j++;
            if (j >= args->max_amp) {
                fprintf(stderr, "[ampliconstats] error: too many amplicons"
                        " (%d). Use -a option to raise this.\n", j);
                return -1;
            }
            amp[j].max_left = 0;
            amp[j].min_right = INT64_MAX;
            amp[j].min_left = INT64_MAX;
            amp[j].max_right = 0;
        }
        if (sites->bp[i].rev == 0) {
            if (i == 0 || last_rev) {
                if (j>0) fprintf(ofp, "\n");
                if (args->multi_ref)
                    fprintf(ofp, "AMPLICON\t%s\t%d", ref, j+1 + first_amp);
                else
                    fprintf(ofp, "AMPLICON\t%d", j+1);
            }
            if (amp[j].nleft >= MAX_PRIMER_PER_AMPLICON) {
                print_error_errno("ampliconstats",
                                  "too many primers per amplicon (%d).\n",
                                  MAX_PRIMER_PER_AMPLICON);
                return -1;
            }
            amp[j].left[amp[j].nleft++] = sites->bp[i].right;
            if (amp[j].max_left < sites->bp[i].right+1)
                amp[j].max_left = sites->bp[i].right+1;
            if (amp[j].min_left > sites->bp[i].right+1)
                amp[j].min_left = sites->bp[i].right+1;
            // BED file, so left+1 as zero based. right(+1-1) as
            // BED goes one beyond end (and we want inclusive range).
            fprintf(ofp, "%c%"PRId64"-%"PRId64, "\t,"[amp[j].nleft > 1],
                    sites->bp[i].left+1, sites->bp[i].right);
        } else {
            if (amp[j].nright >= MAX_PRIMER_PER_AMPLICON) {
                print_error_errno("ampliconstats",
                                  "too many primers per amplicon (%d)",
                                  MAX_PRIMER_PER_AMPLICON);
                return -1;
            }
            amp[j].right[amp[j].nright++] = sites->bp[i].left;
            if (amp[j].min_right > sites->bp[i].left-1)
                amp[j].min_right = sites->bp[i].left-1;
            if (amp[j].max_right < sites->bp[i].left-1) {
                amp[j].max_right = sites->bp[i].left-1;
                if (amp[j].max_right - amp[j].min_left + 1 >=
                    args->max_amp_len) {
                    fprintf(stderr, "[ampliconstats] error: amplicon "
                            "longer (%d) than max_amp_len option (%d)\n",
                            (int)(amp[j].max_right - amp[j].min_left + 2),
                            args->max_amp_len);
                    return -1;
                }
                if (max_right < amp[j].max_right)
                    max_right = amp[j].max_right;
            }
            fprintf(ofp, "%c%"PRId64"-%"PRId64, "\t,"[amp[j].nright > 1],
                    sites->bp[i].left+1, sites->bp[i].right);
        }
        last_rev = sites->bp[i].rev;
    }
    if (last_rev != 1) {
        fprintf(ofp, "\n"); // useful if going to stdout
        fprintf(stderr, "[ampliconstats] error: bed file does not end on"
                " a reverse strand primer.\n");
        return -1;
    }
    *namp = ++j;
    if (j) fprintf(ofp, "\n");

    if (j >= args->max_amp) {
        fprintf(stderr, "[ampliconstats] error: "
                "too many amplicons (%d). Use -a option to raise this.", j);
        return -1;
    }

//    for (i = 0; i < *namp; i++) {
//      printf("%d\t%ld", i, amp[i].length);
//      for (j = 0; j < amp[i].nleft; j++)
//          printf("%c%ld", "\t,"[j>0], amp[i].left[j]);
//      for (j = 0; j < amp[i].nright; j++)
//          printf("%c%ld", "\t,"[j>0], amp[i].right[j]);
//      printf("\n");
//    }

    return max_right;
}

void stats_free(astats_t *st) {
    if (!st)
        return;

    free(st->nreads);
    free(st->nreads2);
    free(st->nfull_reads);
    free(st->nrperc);
    free(st->nrperc2);
    free(st->nbases);
    free(st->nbases2);
    free(st->coverage);
    free(st->covered_perc);
    free(st->covered_perc2);
    free(st->amp_dist);

    free(st->depth_valid);
    free(st->depth_all);

    if (st->tcoord) {
        int i;
        for (i = 0; i <= st->max_amp; i++) {
            if (st->tcoord[i])
                kh_destroy(tcoord, st->tcoord[i]);
        }
        free(st->tcoord);
    }

    khiter_t k;
    for (k = kh_begin(st->qend); k != kh_end(st->qend); k++)
        if (kh_exist(st->qend, k))
            free((void *)kh_key(st->qend, k));
    kh_destroy(qname, st->qend);

    free(st);
}

astats_t *stats_alloc(int64_t max_len, int max_amp, int max_amp_len) {
    astats_t *st = calloc(1, sizeof(*st));
    if (!st)
        return NULL;

    st->max_amp = max_amp;
    st->max_amp_len = max_amp_len;
    st->max_len = max_len;

    if (!(st->nreads  = calloc(max_amp, sizeof(*st->nreads))))  goto err;
    if (!(st->nreads2 = calloc(max_amp, sizeof(*st->nreads2)))) goto err;
    if (!(st->nrperc  = calloc(max_amp, sizeof(*st->nrperc))))  goto err;
    if (!(st->nrperc2 = calloc(max_amp, sizeof(*st->nrperc2)))) goto err;
    if (!(st->nbases  = calloc(max_amp, sizeof(*st->nbases))))  goto err;
    if (!(st->nbases2 = calloc(max_amp, sizeof(*st->nbases2)))) goto err;

    if (!(st->nfull_reads = calloc(max_amp, sizeof(*st->nfull_reads))))
        goto err;

    if (!(st->coverage = calloc(max_amp*max_amp_len, sizeof(*st->coverage))))
        goto err;

    if (!(st->covered_perc  = calloc(max_amp, sizeof(*st->covered_perc))))
        goto err;
    if (!(st->covered_perc2 = calloc(max_amp, sizeof(*st->covered_perc2))))
        goto err;

    if (!(st->tcoord = calloc(max_amp+1, sizeof(*st->tcoord)))) goto err;
    int i;
    for (i = 0; i <= st->max_amp; i++)
        if (!(st->tcoord[i] = kh_init(tcoord)))
            goto err;

    if (!(st->qend = kh_init(qname)))
        goto err;

    if (!(st->depth_valid = calloc(max_len, sizeof(*st->depth_valid))))
        goto err;
    if (!(st->depth_all   = calloc(max_len, sizeof(*st->depth_all))))
        goto err;

    if (!(st->amp_dist  = calloc(max_amp, sizeof(*st->amp_dist))))  goto err;

    return st;

 err:
    stats_free(st);
    return NULL;
}

static void stats_reset(astats_t *st) {
    st->nseq = 0;
    st->nfiltered = 0;
    st->nfailprimer = 0;

    memset(st->nreads,  0, st->max_amp * sizeof(*st->nreads));
    memset(st->nreads2, 0, st->max_amp * sizeof(*st->nreads2));
    memset(st->nfull_reads, 0, st->max_amp * sizeof(*st->nfull_reads));

    memset(st->nrperc,  0, st->max_amp * sizeof(*st->nrperc));
    memset(st->nrperc2, 0, st->max_amp * sizeof(*st->nrperc2));

    memset(st->nbases,  0, st->max_amp * sizeof(*st->nbases));
    memset(st->nbases2, 0, st->max_amp * sizeof(*st->nbases2));

    memset(st->coverage, 0, st->max_amp * st->max_amp_len
           * sizeof(*st->coverage));
    memset(st->covered_perc,  0, st->max_amp * sizeof(*st->covered_perc));
    memset(st->covered_perc2, 0, st->max_amp * sizeof(*st->covered_perc2));

    // Keep the allocated entries as it's likely all files will share
    // the same keys.  Instead we reset counters to zero for common ones
    // and delete rare ones.
    int i;
    for (i = 0; i <= st->max_amp; i++) {
        khiter_t k;
        for (k = kh_begin(st->tcoord[i]);
             k != kh_end(st->tcoord[i]); k++)
            if (kh_exist(st->tcoord[i], k)) {
                if (kh_value(st->tcoord[i], k) < 5)
                    kh_del(tcoord, st->tcoord[i], k);
                else
                    kh_value(st->tcoord[i], k) = 0;
            }
    }

    khiter_t k;
    for (k = kh_begin(st->qend); k != kh_end(st->qend); k++)
        if (kh_exist(st->qend, k))
            free((void *)kh_key(st->qend, k));
    kh_clear(qname, st->qend);

    memset(st->depth_valid, 0, st->max_len * sizeof(*st->depth_valid));
    memset(st->depth_all,   0, st->max_len * sizeof(*st->depth_all));
    memset(st->amp_dist,  0, st->max_amp * sizeof(*st->amp_dist));
}

static void amp_stats_reset(amplicons_t *amps, int nref) {
    int i;
    for (i = 0; i < nref; i++) {
        if (!amps[i].sites)
            continue;
        stats_reset(amps[i].lstats);
    }
}

static int accumulate_stats(astats_args_t *args, amplicons_t *amps,
                            bam1_t *b) {
    int ref = b->core.tid;
    amplicon_t *amp = amps[ref].amp;
    astats_t *stats = amps[ref].lstats;
    int len = amps[ref].len;

    if (!stats)
        return 0;

    stats->nseq++;
    if ((b->core.flag & args->flag_require) != args->flag_require ||
        (b->core.flag & args->flag_filter)  != 0) {
        stats->nfiltered++;
        return 0;
    }

    int64_t start = b->core.pos, mstart = start; // modified start
    int64_t end = bam_endpos(b), i;

    // Compute all-template-depth and valid-template-depth.
    // We track current end location per read name so we can remove overlaps.
    // Potentially we could use this data for a better amplicon-depth
    // count too, but for now it's purely for the per-base plots.
    int ret;
    khiter_t k;
    int prev_start = 0, prev_end = 0;
    if ((b->core.flag & BAM_FPAIRED)
        && !(b->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY))) {
        k = kh_put(qname, stats->qend, bam_get_qname(b), &ret);
        if (ret == 0) {
            prev_start = kh_value(stats->qend, k) & 0xffffffff;
            prev_end = kh_value(stats->qend, k)>>32;
            mstart = MAX(mstart, prev_end);
            // Ideally we'd reuse strings so we don't thrash free/malloc.
            // However let's see if the official way of doing that (malloc
            // itself) is fast enough first.
            free((void *)kh_key(stats->qend, k));
            kh_del(qname, stats->qend, k);
            //fprintf(stderr, "remove overlap %d to %d\n", (int)start, (int)mstart);
        } else {
            if (!(kh_key(stats->qend, k) = strdup(bam_get_qname(b))))
                return -1;

            kh_value(stats->qend, k) = start | (end << 32);
        }
    }
    for (i = mstart; i < end && i < len; i++)
        stats->depth_all[i]++;
    if (i < end) {
        print_error("ampliconstats", "record %s overhangs end of reference",
                    bam_get_qname(b));
        // But keep going, as it's harmless.
    }

    // On single ended runs, eg ONT or PacBio, we just use the start/end
    // of the template to assign.
    int anum = (b->core.flag & BAM_FREVERSE) || !(b->core.flag & BAM_FPAIRED)
        ? (end-1 >= 0 && end-1 < len ? pos2end[end-1] : -1)
        : (start >= 0 && start < len ? pos2start[start] : -1);

    // ivar sometimes soft-clips 100% of the bases.
    // This is essentially unmapped
    if (end == start && (args->flag_filter & BAM_FUNMAP)) {
        stats->nfiltered++;
        return 0;
    }

    if (anum == -1)
        stats->nfailprimer++;

    if (anum >= 0) {
        int64_t c = MIN(end,amp[anum].min_right+1) - MAX(start,amp[anum].max_left);
        if (c > 0) {
            stats->nreads[anum]++;
            // NB: ref bases rather than read bases
            stats->nbases[anum] += c;

            int64_t i;
            if (start < 0) start = 0;
            if (end > len) end = len;

            int64_t ostart = MAX(start, amp[anum].min_left-1);
            int64_t oend = MIN(end, amp[anum].max_right);
            int64_t offset = amp[anum].min_left-1;
            for (i = ostart; i < oend; i++)
                stats->coverage[anum*stats->max_amp_len + i-offset]++;
        } else {
            stats->nfailprimer++;
        }
    }

    // Template length in terms of amplicon number to amplicon number.
    // We expect left to right of same amplicon (len 0), but it may go
    // to next amplicon (len 1) or prev (len -1), etc.
    int64_t t_end;
    int oth_anum = -1;

    if (b->core.flag & BAM_FPAIRED) {
        t_end = (b->core.flag & BAM_FREVERSE ? end : start)
            + b->core.isize;

        // If we've clipped the primers but not followed up with a fixmates
        // then our start+TLEN will take us to a location which is
        // length(LEFT_PRIMER) + length(RIGHT_PRIMER) too far away.
        //
        // The correct solution is to run samtools fixmate so TLEN is correct.
        // The hacky solution is to fudge the expected tlen by double the
        // average primer length (e.g. 50).
        t_end += b->core.isize > 0 ? -args->tlen_adj : +args->tlen_adj;

        if (t_end > 0 && t_end < len && b->core.isize != 0)
            oth_anum = (b->core.flag & BAM_FREVERSE)
                ? pos2start[t_end]
                : pos2end[t_end];
    } else {
        // Not paired (see int anum = (REV || !PAIR) ?en :st expr above)
        oth_anum = pos2start[start];
        t_end = end;
    }

    // We don't want to count our pairs twice.
    // If both left/right are known, count it on left only.
    // If only one is known, we'll only get to this code once
    // so we can also count it.
    int astatus = 2;
    if (anum != -1 && oth_anum != -1) {
        astatus = oth_anum == anum ? 0 : 1;
        if (start <= t_end)
            stats->amp_dist[anum][astatus]++;
    } else if (anum >= 0) {
        stats->amp_dist[anum][astatus = 2]++;
    }

    if (astatus == 0 && !(b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP))) {
        if (prev_end && mstart > prev_end) {
            // 2nd read with gap to 1st; undo previous increment.
            for (i = prev_start; i < prev_end; i++)
                stats->depth_valid[i]--;
            stats->nfull_reads[anum] -= (b->core.flag & BAM_FPAIRED) ? 0.5 : 1;
        } else {
            // 1st read, or 2nd read that overlaps 1st
            for (i = mstart; i < end; i++)
                stats->depth_valid[i]++;
            stats->nfull_reads[anum] += (b->core.flag & BAM_FPAIRED) ? 0.5 : 1;
        }
    }

    // Track template start,end frequencies, so we can give stats on
    // amplicon primer usage.
    if ((b->core.flag & BAM_FPAIRED) && b->core.isize <= 0)
        // left to right only, so we don't double count template positions.
        return 0;

    start = b->core.pos;
    t_end = b->core.flag & BAM_FPAIRED
        ? start + b->core.isize-1
        : end;
    uint64_t tcoord = MIN(start+1, UINT32_MAX) | (MIN(t_end+1, UINT32_MAX)<<32);
    k = kh_put(tcoord, stats->tcoord[anum+1], tcoord, &ret);
    if (ret < 0)
        return -1;
    if (ret == 0)
        kh_value(stats->tcoord[anum+1], k)++;
    else
        kh_value(stats->tcoord[anum+1], k)=1;
    kh_value(stats->tcoord[anum+1], k) |= ((int64_t)astatus<<32);

    return 0;
}

// Append file local stats to global stats
int append_lstats(astats_t *lstats, astats_t *gstats, int namp, int all_nseq) {
    gstats->nseq += lstats->nseq;
    gstats->nfiltered += lstats->nfiltered;
    gstats->nfailprimer += lstats->nfailprimer;

    int a;
    for (a = -1; a < namp; a++) {
        // Add khash local (kl) to khash global (kg)
        khiter_t kl, kg;
        for (kl = kh_begin(lstats->tcoord[a+1]);
             kl != kh_end(lstats->tcoord[a+1]); kl++) {
            if (!kh_exist(lstats->tcoord[a+1], kl) ||
                kh_value(lstats->tcoord[a+1], kl) == 0)
                continue;

            int ret;
            kg = kh_put(tcoord, gstats->tcoord[a+1],
                        kh_key(lstats->tcoord[a+1], kl),
                        &ret);
            if (ret < 0)
                return -1;

            kh_value(gstats->tcoord[a+1], kg) =
                (ret == 0
                 ? (kh_value(gstats->tcoord[a+1], kg) & 0xFFFFFFFF)
                 : 0)
                + kh_value(lstats->tcoord[a+1], kl);
        }
        if (a == -1) continue;

        gstats->nreads[a]  += lstats->nreads[a];
        gstats->nreads2[a] += lstats->nreads[a] * lstats->nreads[a];
        gstats->nfull_reads[a] += lstats->nfull_reads[a];

        // To get mean & sd for amplicon read percentage, we need
        // to do the divisions here as nseq differs for each sample.
        double nrperc = all_nseq ? 100.0 * lstats->nreads[a] / all_nseq : 0;
        gstats->nrperc[a]  += nrperc;
        gstats->nrperc2[a] += nrperc*nrperc;

        gstats->nbases[a]  += lstats->nbases[a];
        gstats->nbases2[a] += lstats->nbases[a] * lstats->nbases[a];

        int d;
        for (d = 0; d < MAX_DEPTH; d++) {
            gstats->covered_perc[a][d]  += lstats->covered_perc[a][d];
            gstats->covered_perc2[a][d] += lstats->covered_perc[a][d]
                                         * lstats->covered_perc[a][d];
        }

        for (d = 0; d < 3; d++)
            gstats->amp_dist[a][d] += lstats->amp_dist[a][d];
    }

    for (a = 0; a < lstats->max_len; a++) {
        gstats->depth_valid[a] += lstats->depth_valid[a];
        gstats->depth_all[a]   += lstats->depth_all[a];
    }

    return 0;
}

int append_stats(amplicons_t *amps, int nref) {
    int i, r, all_nseq = 0;
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = amps[r].lstats;
        all_nseq  += stats->nseq - stats->nfiltered - stats->nfailprimer;
    }

    for (i = 0; i < nref; i++) {
        if (!amps[i].sites)
            continue;
        if (append_lstats(amps[i].lstats, amps[i].gstats, amps[i].namp,
                          all_nseq) < 0)
            return -1;
    }

    return 0;
}

typedef struct {
    int32_t start, end;
    uint32_t freq;
    uint32_t status;
} tcoord_t;

// Sort tcoord by descending frequency and then ascending start and  end.
static int tcoord_freq_sort(const void *vp1, const void *vp2) {
    const tcoord_t *t1 = (const tcoord_t *)vp1;
    const tcoord_t *t2 = (const tcoord_t *)vp2;

    if (t1->freq != t2->freq)
        return t2->freq - t1->freq;

    if (t1->start != t2->start)
        return t1->start - t2->start;

    return t1->end - t2->end;
}


/*
 * Merges tcoord start,end,freq,status tuples if their coordinates are
 * close together.  We aim to keep the start,end for the most frequent
 * value and assume that is the correct coordinate and all others are
 * minor fluctuations due to errors or variants.
 *
 * We sort by frequency first and then merge later items in the list into
 * the earlier more frequent ones.  It's O(N^2), but sufficient for now
 * given current scale of projects.
 *
 * If we ever need to resolve that then consider sorting by start
 * coordinate and scanning the list to find all items within X, find
 * the most frequent of those, and then cluster that way.  (I'd have
 * done that had I thought of it at the time!)
 */
static void aggregate_tcoord(astats_args_t *args, tcoord_t *tpos, size_t *np){
    size_t n = *np, j, j2, j3, k;

    // Sort by frequency and cluster infrequent coords into frequent
    // ones provided they're close by.
    // This is O(N^2), but we've already binned by tcoord_bin/2 so
    // the list isn't intended to be vast at this point.
    qsort(tpos, n, sizeof(*tpos), tcoord_freq_sort);

    // For frequency ties, find mid start coord, and then find mid end
    // coord of those matching start.
    // We make that the first item so we merge into that mid point.
    for (j = 0; j < n; j++) {
        for (j2 = j+1; j2 < n; j2++) {
            if (tpos[j].freq != tpos[j2].freq)
                break;
            if (tpos[j2].start - tpos[j].start >= args->tcoord_bin)
                break;
        }

        // j to j2 all within bin of a common start,
        // m is the mid start.
        if (j2-1 > j) {
            size_t m = (j2-1 + j)/2;

            // Find mid end for this same start
            while (m > 1 && tpos[m].start == tpos[m-1].start)
                m--;
            for (j3 = m+1; j3 < j2; j3++) {
                if (tpos[m].start != tpos[j3].start)
                    break;
                if (tpos[m].end - tpos[j3].end >= args->tcoord_bin)
                    break;
            }
            if (j3-1 > m)
                m = (j3-1 + m)/2;

            // Swap with first item.
            tcoord_t tmp = tpos[j];
            tpos[j] = tpos[m];
            tpos[m] = tmp;
            j = j2-1;
        }
    }

    // Now merge in coordinates.
    // This bit is O(N^2), so consider binning first to reduce the
    // size of the list if we have excessive positional variation.
    for (k = j = 0; j < n; j++) {
        if (!tpos[j].freq)
            continue;

        if (k < j)
            tpos[k] = tpos[j];

        for (j2 = j+1; j2 < n; j2++) {
            if (ABS(tpos[j].start-tpos[j2].start) < args->tcoord_bin/2 &&
                ABS(tpos[j].end  -tpos[j2].end)  < args->tcoord_bin/2 &&
                tpos[j].status == tpos[j2].status) {
                tpos[k].freq += tpos[j2].freq;
                tpos[j2].freq = 0;
            }
        }
        k++;
    }

    *np = k;
}

int dump_stats(astats_args_t *args, char type, char *name, int nfile,
               amplicons_t *amps, int nref, int local) {
    int i, r;
    FILE *ofp = args->out_fp;
    tcoord_t *tpos = NULL;
    size_t ntcoord = 0;

    // summary stats for this sample (or for all samples)
    fprintf(ofp, "# Summary stats.\n");
    fprintf(ofp, "# Use 'grep ^%cSS | cut -f 2-' to extract this part.\n", type);

    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        int nmatch = stats->nseq - stats->nfiltered - stats->nfailprimer;
        char *name_ref = malloc(strlen(name) + strlen(amps[r].ref) + 2);
        if (!name_ref)
            return -1;
        if (args->multi_ref)
            sprintf(name_ref, "%s\t%s", name, amps[r].ref);
        else
            sprintf(name_ref, "%s", name);
        fprintf(ofp, "%cSS\t%s\traw total sequences:\t%d\n",
                type, name_ref, stats->nseq);
        fprintf(ofp, "%cSS\t%s\tfiltered sequences:\t%d\n",
                type, name_ref, stats->nfiltered);
        fprintf(ofp, "%cSS\t%s\tfailed primer match:\t%d\n",
                type, name_ref, stats->nfailprimer);
        fprintf(ofp, "%cSS\t%s\tmatching sequences:\t%d\n",
                type, name_ref, nmatch);

        int d = 0;
        do {
            // From first to last amplicon only, so not entire consensus.
            // If contig length is known, maybe we want to add the missing
            // count to < DEPTH figures?
            int64_t start = 0, covered = 0, total = 0;
            amplicon_t *amp = amps[r].amp;
            for (i = 0; i < amps[r].namp; i++) {
                int64_t j, offset = amp[i].min_left-1;
                if (amp[i].min_right - amp[i].min_left > stats->max_amp_len) {
                    fprintf(stderr, "[ampliconstats] error: "
                            "Maximum amplicon length (%d) exceeded for '%s'\n",
                            stats->max_amp, name);
                    return -1;
                }
                for (j = MAX(start, amp[i].max_left-1);
                     j < MAX(start, amp[i].min_right); j++) {
                    if (stats->coverage[i*stats->max_amp_len + j-offset]
                        >= args->min_depth[d])
                        covered++;
                    total++;
                }
                start = MAX(start, amp[i].min_right);
            }
            fprintf(ofp, "%cSS\t%s\tconsensus depth count < %d and >= %d:\t%"
                    PRId64"\t%"PRId64"\n", type, name_ref,
                    args->min_depth[d], args->min_depth[d],
                    total-covered, covered);
        } while (++d < MAX_DEPTH && args->min_depth[d]);

        free(name_ref);
    }

    // Read count
    fprintf(ofp, "# Absolute matching read counts per amplicon.\n");
    fprintf(ofp, "# Use 'grep ^%cREADS | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cREADS\t%s", type, name);
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        for (i = 0; i < amps[r].namp; i++) {
            fprintf(ofp, "\t%"PRId64, stats->nreads[i]);
        }
    }
    fprintf(ofp, "\n");

    // Valid depth is the number of full length reads (already divided
    // by the number we expect to cover), so +0.5 per read in pair.
    // A.k.a "usable depth" in the plots.
    fprintf(ofp, "%cVDEPTH\t%s", type, name);
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        for (i = 0; i < amps[r].namp; i++)
            fprintf(ofp, "\t%d", (int)stats->nfull_reads[i]);
    }
    fprintf(ofp, "\n");

    if (type == 'C') {
        // For combined we can compute mean & standard deviation too
        fprintf(ofp, "CREADS\tMEAN");
        for (r = 0; r < nref; r++) {
            if (!amps[r].sites)
                continue;
            astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
            for (i = 0; i < amps[r].namp; i++) {
                fprintf(ofp, "\t%.1f", stats->nreads[i] / (double)nfile);
            }
        }
        fprintf(ofp, "\n");

        fprintf(ofp, "CREADS\tSTDDEV");
        for (r = 0; r < nref; r++) {
            if (!amps[r].sites)
                continue;
            astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
            for (i = 0; i < amps[r].namp; i++) {
                double n1 = stats->nreads[i];
                fprintf(ofp, "\t%.1f", nfile > 1 && stats->nreads2[i] > 0
                        ? sqrt(stats->nreads2[i]/(double)nfile
                               - (n1/nfile)*(n1/nfile))
                        : 0);
            }
        }
        fprintf(ofp, "\n");
    }

    fprintf(ofp, "# Read percentage of distribution between amplicons.\n");
    fprintf(ofp, "# Use 'grep ^%cRPERC | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cRPERC\t%s", type, name);
    int all_nseq = 0;
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        all_nseq  += stats->nseq - stats->nfiltered - stats->nfailprimer;
    }
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        for (i = 0; i < amps[r].namp; i++) {
            if (type == 'C') {
                fprintf(ofp, "\t%.3f", (double)stats->nrperc[i] / nfile);
            } else {
                fprintf(ofp, "\t%.3f",
                        all_nseq ? 100.0 * stats->nreads[i] / all_nseq : 0);
            }
        }
    }
    fprintf(ofp, "\n");

    if (type == 'C') {
        // For combined we compute mean and standard deviation too
        fprintf(ofp, "CRPERC\tMEAN");
        for (r = 0; r < nref; r++) {
            if (!amps[r].sites)
                continue;
            astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
            for (i = 0; i < amps[r].namp; i++) {
                fprintf(ofp, "\t%.3f", stats->nrperc[i] / nfile);
            }
        }
        fprintf(ofp, "\n");

        fprintf(ofp, "CRPERC\tSTDDEV");
        for (r = 0; r < nref; r++) {
            if (!amps[r].sites)
                continue;
            astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
            for (i = 0; i < amps[r].namp; i++) {
                // variance = SUM(X^2) - ((SUM(X)^2) / N)
                double n1 = stats->nrperc[i];
                double v = stats->nrperc2[i]/nfile - (n1/nfile)*(n1/nfile);
                fprintf(ofp, "\t%.3f", v>0?sqrt(v):0);
            }
        }
        fprintf(ofp, "\n");
    }

    // Base depth
    fprintf(ofp, "# Read depth per amplicon.\n");
    fprintf(ofp, "# Use 'grep ^%cDEPTH | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cDEPTH\t%s", type, name);
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        amplicon_t *amp = amps[r].amp;
        for (i = 0; i < amps[r].namp; i++) {
            int nseq = stats->nseq - stats->nfiltered - stats->nfailprimer;
            int64_t alen = amp[i].min_right - amp[i].max_left+1;
            fprintf(ofp, "\t%.1f", nseq ? stats->nbases[i] / (double)alen : 0);
        }
    }
    fprintf(ofp, "\n");

    if (type == 'C') {
        // For combined we can compute mean & standard deviation too
        fprintf(ofp, "CDEPTH\tMEAN");
        for (r = 0; r < nref; r++) {
            if (!amps[r].sites)
                continue;
            astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
            amplicon_t *amp = amps[r].amp;
            int nseq = stats->nseq - stats->nfiltered - stats->nfailprimer;
            for (i = 0; i < amps[r].namp; i++) {
                int64_t alen = amp[i].min_right - amp[i].max_left+1;
                fprintf(ofp, "\t%.1f", nseq ? stats->nbases[i] / (double)alen / nfile : 0);
            }
        }
        fprintf(ofp, "\n");

        fprintf(ofp, "CDEPTH\tSTDDEV");
        for (r = 0; r < nref; r++) {
            if (!amps[r].sites)
                continue;
            astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
            amplicon_t *amp = amps[r].amp;
            for (i = 0; i < amps[r].namp; i++) {
                double alen = amp[i].min_right - amp[i].max_left+1;
                double n1 = stats->nbases[i] / alen;
                double v = stats->nbases2[i] / (alen*alen) /nfile
                    - (n1/nfile)*(n1/nfile);
                fprintf(ofp, "\t%.1f", v>0?sqrt(v):0);
            }
        }
        fprintf(ofp, "\n");
    }

    // Percent Coverage
    if (type == 'F') {
        fprintf(ofp, "# Percentage coverage per amplicon\n");
        fprintf(ofp, "# Use 'grep ^%cPCOV | cut -f 2-' to extract this part.\n", type);
        int d = 0;
        do {
            fprintf(ofp, "%cPCOV-%d\t%s", type, args->min_depth[d], name);

            for (r = 0; r < nref; r++) {
                if (!amps[r].sites)
                    continue;
                astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
                amplicon_t *amp = amps[r].amp;
                for (i = 0; i < amps[r].namp; i++) {
                    int covered = 0;
                    if (amp[i].min_right - amp[i].min_left > stats->max_amp_len) {
                        fprintf(stderr, "[ampliconstats] error: "
                                "Maximum amplicon length (%d) exceeded for '%s'\n",
                                stats->max_amp, name);
                        return -1;
                    }
                    int64_t j, offset = amp[i].min_left-1;
                    for (j = amp[i].max_left-1; j < amp[i].min_right; j++) {
                        int apos = i*stats->max_amp_len + j-offset;
                        if (stats->coverage[apos] >= args->min_depth[d])
                            covered++;
                    }
                    int64_t alen = amp[i].min_right - amp[i].max_left+1;
                    stats->covered_perc[i][d] = 100.0 * covered / alen;
                    fprintf(ofp, "\t%.2f", 100.0 * covered / alen);
                }
            }
            fprintf(ofp, "\n");
        } while (++d < MAX_DEPTH && args->min_depth[d]);

    } else if (type == 'C') {
        // For combined we can compute mean & standard deviation too
        int d = 0;
        do {
            fprintf(ofp, "CPCOV-%d\tMEAN", args->min_depth[d]);
            for (r = 0; r < nref; r++) {
                if (!amps[r].sites)
                    continue;
                astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
                for (i = 0; i < amps[r].namp; i++) {
                    fprintf(ofp, "\t%.1f", stats->covered_perc[i][d] / nfile);
                }
            }
            fprintf(ofp, "\n");

            fprintf(ofp, "CPCOV-%d\tSTDDEV", args->min_depth[d]);
            for (r = 0; r < nref; r++) {
                if (!amps[r].sites)
                    continue;
                astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
                for (i = 0; i < amps[r].namp; i++) {
                    double n1 = stats->covered_perc[i][d] / nfile;
                    double v = stats->covered_perc2[i][d] / nfile - n1*n1;
                    fprintf(ofp, "\t%.1f", v>0?sqrt(v):0);
                }
            }
            fprintf(ofp, "\n");
        } while (++d < MAX_DEPTH && args->min_depth[d]);
    }

    // Plus base depth for all reads, irrespective of amplicon.
    // This is post overlap removal, if reads in the read-pair overlap.
    fprintf(ofp, "# Depth per reference base for ALL data.\n");
    fprintf(ofp, "# Use 'grep ^%cDP_ALL | cut -f 2-' to extract this part.\n",
            type);
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        if (args->multi_ref)
            fprintf(ofp, "%cDP_ALL\t%s\t%s", type, name, amps[r].ref);
        else
            fprintf(ofp, "%cDP_ALL\t%s", type, name);

        for (i = 0; i < amps[r].len; i++) {
            // Basic run-length encoding provided all values are within
            // +- depth_bin fraction of the mid-point.
            int dmin = stats->depth_all[i], dmax = stats->depth_all[i], j;
            double dmid = (dmin + dmax)/2.0;
            double low  = dmid*(1-args->depth_bin);
            double high = dmid*(1+args->depth_bin);
            for (j = i+1; j < amps[r].len; j++) {
                int d = stats->depth_all[j];
                if (d < low || d > high)
                    break;
                if (dmin > d) {
                    dmin = d;
                    dmid = (dmin + dmax)/2.0;
                    low  = dmid*(1-args->depth_bin);
                    high = dmid*(1+args->depth_bin);
                } else if (dmax < d) {
                    dmax = d;
                    dmid = (dmin + dmax)/2.0;
                    low  = dmid*(1-args->depth_bin);
                    high = dmid*(1+args->depth_bin);
                }
            }
            fprintf(ofp, "\t%d,%d", (int)dmid, j-i);
            i = j-1;
        }
        fprintf(ofp, "\n");
    }

    // And depth for only reads matching to a single amplicon for full
    // length.  This is post read overlap removal.
    fprintf(ofp, "# Depth per reference base for full-length valid amplicon data.\n");
    fprintf(ofp, "# Use 'grep ^%cDP_VALID | cut -f 2-' to extract this "
            "part.\n", type);
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        if (args->multi_ref)
            fprintf(ofp, "%cDP_VALID\t%s\t%s", type, name, amps[r].ref);
        else
            fprintf(ofp, "%cDP_VALID\t%s", type, name);

        for (i = 0; i < amps[r].len; i++) {
            int dmin = stats->depth_valid[i], dmax = stats->depth_valid[i], j;
            double dmid = (dmin + dmax)/2.0;
            double low  = dmid*(1-args->depth_bin);
            double high = dmid*(1+args->depth_bin);
            for (j = i+1; j < amps[r].len; j++) {
                int d = stats->depth_valid[j];
                if (d < low || d > high)
                    break;
                if (dmin > d) {
                    dmin = d;
                    dmid = (dmin + dmax)/2.0;
                    low  = dmid*(1-args->depth_bin);
                    high = dmid*(1+args->depth_bin);
                } else if (dmax < d) {
                    dmax = d;
                    dmid = (dmin + dmax)/2.0;
                    low  = dmid*(1-args->depth_bin);
                    high = dmid*(1+args->depth_bin);
                }
            }
            fprintf(ofp, "\t%d,%d", (int)dmid, j-i);
            i = j-1;
        }
        fprintf(ofp, "\n");
    }

    // TCOORD (start to end) distribution
    fprintf(ofp, "# Distribution of aligned template coordinates.\n");
    fprintf(ofp, "# Use 'grep ^%cTCOORD | cut -f 2-' to extract this part.\n", type);
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        for (i = 0 - (nref==1); i < amps[r].namp; i++) {
            if (ntcoord < kh_size(stats->tcoord[i+1])) {
                ntcoord = kh_size(stats->tcoord[i+1]);
                tcoord_t *tmp = realloc(tpos, ntcoord * sizeof(*tmp));
                if (!tmp) {
                    free(tpos);
                    return -1;
                }
                tpos = tmp;
            }

            khiter_t k;
            size_t n = 0, j;
            for (k = kh_begin(stats->tcoord[i+1]);
                 k != kh_end(stats->tcoord[i+1]); k++) {
                if (!kh_exist(stats->tcoord[i+1], k) ||
                    (kh_value(stats->tcoord[i+1], k) & 0xFFFFFFFF) == 0)
                    continue;
                // Key is start,end in 32-bit quantities.
                // Yes this limits us to 4Gb references, but just how
                // many primers are we planning on making?  Not that many
                // I hope.
                tpos[n].start = kh_key(stats->tcoord[i+1], k)&0xffffffff;
                tpos[n].end   = kh_key(stats->tcoord[i+1], k)>>32;

                // Value is frequency (top 32-bits) and status (bottom 32).
                tpos[n].freq   = kh_value(stats->tcoord[i+1], k)&0xffffffff;
                tpos[n].status = kh_value(stats->tcoord[i+1], k)>>32;
                n++;
            }

            if (args->tcoord_bin > 1)
                aggregate_tcoord(args, tpos, &n);

            fprintf(ofp, "%cTCOORD\t%s\t%d", type, name,
                    i+1+amps[r].first_amp); // per amplicon
            for (j = 0; j < n; j++) {
                if (tpos[j].freq < args->tcoord_min_count)
                    continue;
                fprintf(ofp, "\t%d,%d,%u,%u",
                        tpos[j].start,
                        tpos[j].end,
                        tpos[j].freq,
                        tpos[j].status);
            }
            fprintf(ofp, "\n");
        }
    }


    // AMP length distribution.
    // 0 = both ends in this amplicon
    // 1 = ends in different amplicons
    // 2 = other end matching an unknown amplicon site
    //     (see tcoord for further analysis of where)
    fprintf(ofp, "# Classification of amplicon status.  Columns are\n");
    fprintf(ofp, "# number with both primers from this amplicon, number with\n");
    fprintf(ofp, "# primers from different amplicon, and number with a position\n");
    fprintf(ofp, "# not matching any valid amplicon primer site\n");
    fprintf(ofp, "# Use 'grep ^%cAMP | cut -f 2-' to extract this part.\n", type);

    fprintf(ofp, "%cAMP\t%s\t0", type, name); // all merged
    int amp_dist[3] = {0};
    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        for (i = 0; i < amps[r].namp; i++) { // accumulate for all amps
            amp_dist[0] += stats->amp_dist[i][0];
            amp_dist[1] += stats->amp_dist[i][1];
            amp_dist[2] += stats->amp_dist[i][2];
        }
    }
    fprintf(ofp, "\t%d\t%d\t%d\n", amp_dist[0], amp_dist[1], amp_dist[2]);

    for (r = 0; r < nref; r++) {
        if (!amps[r].sites)
            continue;
        astats_t *stats = local ? amps[r].lstats : amps[r].gstats;
        for (i = 0; i < amps[r].namp; i++) {
            // per amplicon
            fprintf(ofp, "%cAMP\t%s\t%d", type, name, i+1+amps[r].first_amp);
            fprintf(ofp, "\t%d\t%d\t%d\n", stats->amp_dist[i][0],
                    stats->amp_dist[i][1], stats->amp_dist[i][2]);
        }
    }

    free(tpos);
    return 0;
}

int dump_lstats(astats_args_t *args, char type, char *name, int nfile,
               amplicons_t *amps, int nref) {
    return dump_stats(args, type, name, nfile, amps, nref, 1);
}

int dump_gstats(astats_args_t *args, char type, char *name, int nfile,
               amplicons_t *amps, int nref) {
    return dump_stats(args, type, name, nfile, amps, nref, 0);
}

char const *get_sample_name(sam_hdr_t *header, char *RG) {
    kstring_t ks = {0};
    sam_hdr_find_tag_id(header, "RG", RG?"ID":NULL, RG, "SM", &ks);
    return ks.s;
}

// Return maximum reference length (SQ is NULL) or the length
// of the specified reference in SQ.
int64_t get_ref_len(sam_hdr_t *header, const char *SQ) {
    if (SQ) {
        int tid = SQ ? sam_hdr_name2tid(header, SQ) : 0;
        return tid >= 0 ? sam_hdr_tid2len(header, tid) : -1;
    } else {
        int nref = sam_hdr_nref(header), tid;;
        int64_t len = 0;
        for (tid = 0; tid < nref; tid++) {
            int64_t rl = sam_hdr_tid2len(header, tid);
            if (len < rl)
                len = rl;
        }
        return len;
    }
}

static int amplicon_stats(astats_args_t *args,
                          khash_t(bed_list_hash) *bed_hash,
                          char **filev, int filec) {
    int i, ref = -1, ref_tid = -1, ret = -1, nref = 0;
    samFile *fp = NULL;
    sam_hdr_t *header = NULL;
    bam1_t *b = bam_init1();
    FILE *ofp = args->out_fp;
    char sname_[8192], *sname = NULL;
    amplicons_t *amps = NULL;

    // Report initial SS header.  We gather data from the bed_hash entries
    // as well as from the first SAM header (with the requirement that all
    // headers should be compatible).
    if (filec) {
        if (!(fp = sam_open_format(filev[0], "r", &args->ga.in))) {
            print_error_errno("ampliconstats",
                              "Cannot open input file \"%s\"",
                              filev[0]);
            goto err;
        }
        if (!(header = sam_hdr_read(fp)))
            goto err;

        if (!amps) {
            amps = calloc(nref=sam_hdr_nref(header), sizeof(*amps));
            if (!amps)
                goto err;
            fprintf(ofp, "# Summary statistics, used for scaling the plots.\n");
            fprintf(ofp, "SS\tSamtools version: %s\n", samtools_version());
            fprintf(ofp, "SS\tCommand line: %s\n", args->argv);
            fprintf(ofp, "SS\tNumber of files:\t%d\n", filec);

            // Note: order of hash entries will be different to order of
            // BED file which may also differ to order of SQ headers.
            // SQ header is canonical ordering (pos sorted file).
            khiter_t k;
            int bam_nref = sam_hdr_nref(header);
            for (i = 0; i < bam_nref; i++) {
                k = kh_get(bed_list_hash, bed_hash,
                           sam_hdr_tid2name(header, i));
                if (!kh_exist(bed_hash, k))
                    continue;

                bed_entry_list_t *sites = &kh_value(bed_hash, k);

                ref = i;
                amps[ref].ref = kh_key(bed_hash, k);
                amps[ref].sites = sites;
                amps[ref].namp = count_amplicon(sites);
                amps[ref].amp  = calloc(sites->length,
                                        sizeof(*amps[ref].amp));
                if (!amps[ref].amp)
                    goto err;
                if (args->multi_ref)
                    fprintf(ofp, "SS\tNumber of amplicons:\t%s\t%d\n",
                            kh_key(bed_hash, k), amps[ref].namp);
                else
                    fprintf(ofp, "SS\tNumber of amplicons:\t%d\n",
                            amps[ref].namp);

                amps[ref].tid = ref;
                if (ref_tid == -1)
                    ref_tid = ref;

                int64_t len = get_ref_len(header, kh_key(bed_hash, k));
                amps[ref].len = len;
                if (args->multi_ref)
                    fprintf(ofp, "SS\tReference length:\t%s\t%"PRId64"\n",
                            kh_key(bed_hash, k), len);
                else
                    fprintf(ofp, "SS\tReference length:\t%"PRId64"\n",
                            len);

                amps[ref].lstats = stats_alloc(len, args->max_amp,
                                               args->max_amp_len);
                amps[ref].gstats = stats_alloc(len, args->max_amp,
                                               args->max_amp_len);
                if (!amps[ref].lstats || !amps[ref].gstats)
                    goto err;
            }
        }

        sam_hdr_destroy(header);
        header = NULL;
        if (sam_close(fp) < 0) {
            fp = NULL;
            goto err;
        }
        fp = NULL;
    }
    fprintf(ofp, "SS\tEnd of summary\n");

    // Extract the bits of amplicon data we need from bed hash and turn
    // it into a position-to-amplicon lookup table.
    int offset = 0;
    for (i = 0; i < nref; i++) {
        if (!amps[i].sites)
            continue;

        amps[i].first_amp = offset;
        if (bed2amplicon(args, amps[i].sites, amps[i].amp,
                         &amps[i].namp, i==0, amps[i].ref, offset) < 0)
            goto err;

        offset += amps[i].namp; // cumulative amplicon number across refs
    }

    // Now iterate over file contents, one at a time.
    for (i = 0; i < filec; i++) {
        char *nstart = filev[i];

        fp = sam_open_format(filev[i], "r", &args->ga.in);
        if (!fp) {
            print_error_errno("ampliconstats",
                              "Cannot open input file \"%s\"",
                              filev[i]);
            goto err;
        }

        if (args->ga.nthreads > 0)
            hts_set_threads(fp, args->ga.nthreads);

        if (!(header = sam_hdr_read(fp)))
            goto err;

        if (nref != sam_hdr_nref(header)) {
            print_error_errno("ampliconstats",
                              "SAM headers are not consistent across input files");
            goto err;
        }
        int r;
        for (r = 0; r < nref; r++) {
            if (!amps[r].sites)
                continue;
            if (!amps[r].ref ||
                strcmp(amps[r].ref, sam_hdr_tid2name(header, r)) != 0 ||
                amps[r].len != sam_hdr_tid2len(header, r)) {
                print_error_errno("ampliconstats",
                                  "SAM headers are not consistent across "
                                  "input files");
                goto err;
            }
        }

        if (args->use_sample_name)
            sname = (char *)get_sample_name(header, NULL);

        if (!sname) {
            sname = sname_;
            char *nend = filev[i] + strlen(filev[i]), *cp;
            if ((cp = strrchr(filev[i], '/')))
                nstart = cp+1;
            if ((cp = strrchr(nstart, '.')) &&
                (strcmp(cp, ".bam") == 0 ||
                 strcmp(cp, ".sam") == 0 ||
                 strcmp(cp, ".cram") == 0))
                nend = cp;
            if (nend - nstart >= 8192) nend = nstart+8191;
            memcpy(sname, nstart, nend-nstart);
            sname[nend-nstart] = 0;
        }

        // Stats local to this sample only
        amp_stats_reset(amps, nref);

        int last_ref = -9;
        while ((r = sam_read1(fp, header, b)) >= 0) {
            // Other filter options useful here?
            if (b->core.tid < 0)
                continue;

            if (last_ref != b->core.tid) {
                last_ref  = b->core.tid;
                if (initialise_amp_pos_lookup(args, amps, last_ref) < 0)
                    goto err;
            }

            if (accumulate_stats(args, amps, b) < 0)
                goto err;
        }

        if (r < -1) {
            print_error_errno("ampliconstats", "Fail reading record");
            goto err;
        }

        sam_hdr_destroy(header);
        if (sam_close(fp) < 0) {
            fp = NULL;
            goto err;
        }

        fp = NULL;
        header = NULL;

        if (dump_lstats(args, 'F', sname, filec, amps, nref) < 0)
            goto err;

        if (append_stats(amps, nref) < 0)
            goto err;

        if (sname && sname != sname_)
            free(sname);
        sname = NULL;
    }

    if (dump_gstats(args, 'C', "COMBINED", filec, amps, nref) < 0)
        goto err;

    ret = 0;
 err:
    bam_destroy1(b);
    if (ret) {
        if (header)
            sam_hdr_destroy(header);
        if (fp)
            sam_close(fp);
    }
    for (i = 0; i < nref; i++) {
        stats_free(amps[i].lstats);
        stats_free(amps[i].gstats);
        free(amps[i].amp);
    }
    free(amps);
    free(pos2start);
    free(pos2end);
    if (ret) {
        if (sname && sname != sname_)
            free(sname);
    }

    return ret;
}

static int usage(astats_args_t *args, FILE *fp, int exit_status) {
    fprintf(fp,
"\n"
"Usage: samtools ampliconstats [options] primers.bed *.bam > astats.txt\n"
"\n"
"Options:\n");
    fprintf(fp, "  -f, --required-flag STR|INT\n"
            "               Only include reads with all of the FLAGs present [0x%X]\n",args->flag_require);
    fprintf(fp, "  -F, --filter-flag STR|INT\n"
            "               Only include reads with none of the FLAGs present [0x%X]\n",args->flag_filter & 0xffff);
    fprintf(fp, "  -a, --max-amplicons INT\n"
            "               Change the maximum number of amplicons permitted [%d]\n", MAX_AMP);
    fprintf(fp, "  -l, --max-amplicon-length INT\n"
            "               Change the maximum length of an individual amplicon [%d]\n", MAX_AMP_LEN);
    fprintf(fp, "  -d, --min-depth INT[,INT]...\n"
            "               Minimum base depth(s) to consider position covered [%d]\n", args->min_depth[0]);
    fprintf(fp, "  -m, --pos-margin INT\n"
            "               Margin of error for matching primer positions [%d]\n", args->max_delta);
    fprintf(fp, "  -o, --output FILE\n"
            "               Specify output file [stdout if unset]\n");
    fprintf(fp, "  -s, --use-sample-name\n"
            "               Use the sample name from the first @RG header line\n");
    fprintf(fp, "  -t, --tlen-adjust INT\n"
            "               Add/subtract from TLEN; use when clipping but no fixmate step\n");
    fprintf(fp, "  -b, --tcoord-bin INT\n"
            "               Bin template start,end positions into multiples of INT[1]\n");
    fprintf(fp, "  -c, --tcoord-min-count INT\n"
            "               Minimum template start,end frequency for recording [%d]\n", TCOORD_MIN_COUNT);
    fprintf(fp, "  -D, --depth-bin FRACTION\n"
            "               Merge FDP values within +/- FRACTION together\n");
    fprintf(fp, "  -S, --single-ref\n"
            "               Force single-ref (<=1.12) output format\n");
    sam_global_opt_help(fp, "I.--.@");

    return exit_status;
}

int main_ampliconstats(int argc, char **argv) {
    astats_args_t args = {
        .ga = SAM_GLOBAL_ARGS_INIT,
        .flag_require = 0,
        .flag_filter = 0x10B04,
        //.sites = BED_LIST_INIT,
        .max_delta = 30, // large enough to cope with alt primers
        .min_depth = {1},
        .use_sample_name = 0,
        .max_amp = MAX_AMP,
        .max_amp_len = MAX_AMP_LEN,
        .tlen_adj = 0,
        .out_fp = stdout,
        .tcoord_min_count = TCOORD_MIN_COUNT,
        .tcoord_bin = 1,
        .depth_bin = 0.01,
        .multi_ref = 1
    }, oargs = args;

    static const struct option loptions[] =
    {
        SAM_OPT_GLOBAL_OPTIONS('I', 0, '-', '-', 0, '@'),
        {"help", no_argument, NULL, 'h'},
        {"flag-require", required_argument, NULL, 'f'},
        {"flag-filter", required_argument, NULL, 'F'},
        {"min-depth", required_argument, NULL, 'd'},
        {"output", required_argument, NULL, 'o'},
        {"pos-margin", required_argument, NULL, 'm'},
        {"use-sample-name", no_argument, NULL, 's'},
        {"max-amplicons", required_argument, NULL, 'a'},
        {"max-amplicon-length", required_argument, NULL, 'l'},
        {"tlen-adjust", required_argument, NULL, 't'},
        {"tcoord-min-count", required_argument, NULL, 'c'},
        {"tcoord-bin", required_argument, NULL, 'b'},
        {"depth-bin", required_argument, NULL, 'D'},
        {"single-ref", no_argument, NULL, 'S'},
        {NULL, 0, NULL, 0}
    };
    int opt;

    while ( (opt=getopt_long(argc,argv,"?hf:F:@:p:m:d:sa:l:t:o:c:b:D:S",loptions,NULL))>0 ) {
        switch (opt) {
        case 'f': args.flag_require = bam_str2flag(optarg); break;
        case 'F':
            if (args.flag_filter & 0x10000)
                args.flag_filter = 0; // strip default on first -F usage
            args.flag_filter |= bam_str2flag(optarg); break;

        case 'm': args.max_delta = atoi(optarg); break; // margin
        case 'D': args.depth_bin = atof(optarg); break; // depth bin fraction
        case 'd': {
            int d = 0;
            char *cp = optarg, *ep;
            do {
                long n = strtol(cp, &ep, 10);
                args.min_depth[d++] = n;
                if (*ep != ',')
                    break;
                cp = ep+1;
            } while (d < MAX_DEPTH);
            break;
        }

        case 'a': args.max_amp = atoi(optarg)+1;break;
        case 'l': args.max_amp_len = atoi(optarg)+1;break;

        case 'c': args.tcoord_min_count = atoi(optarg);break;
        case 'b':
            args.tcoord_bin = atoi(optarg);
            if (args.tcoord_bin < 1)
                args.tcoord_bin = 1;
            break;

        case 't': args.tlen_adj = atoi(optarg);break;

        case 's': args.use_sample_name = 1;break;

        case 'o':
            if (!(args.out_fp = fopen(optarg, "w"))) {
                perror(optarg);
                return 1;
            }
            break;

        case 'S':
            args.multi_ref = 0;
            break;

        case '?': return usage(&oargs, stderr, EXIT_FAILURE);
        case 'h': return usage(&oargs, stdout, EXIT_SUCCESS);

        default:
            if (parse_sam_global_opt(opt, optarg, loptions, &args.ga) != 0)
                usage(&oargs,stderr, EXIT_FAILURE);
            break;
        }
    }

    if (argc <= optind)
        return usage(&oargs, stdout, EXIT_SUCCESS);
    if (argc <= optind+1 && isatty(STDIN_FILENO))
        return usage(&oargs, stderr, EXIT_FAILURE);

    khash_t(bed_list_hash) *bed_hash = kh_init(bed_list_hash);
    if (load_bed_file_multi_ref(argv[optind], 1, 0, bed_hash)) {
        print_error_errno("ampliconstats",
                          "Could not read file \"%s\"", argv[optind]);
        return 1;

    }

    khiter_t k, ref_count = 0;
    for (k = kh_begin(bed_hash); k != kh_end(bed_hash); k++) {
        if (!kh_exist(bed_hash, k))
            continue;
        ref_count++;
    }
    if (ref_count == 0)
        return 1;
    if (ref_count > 1 && args.multi_ref == 0) {
        print_error("ampliconstats",
                    "Single-ref mode is not permitted for BED files\n"
                    "containing more than one reference.");
        return 1;
    }

    args.argv = stringify_argv(argc, argv);
    int ret;
    if (argc == ++optind) {
        char *av = "-";
        ret = amplicon_stats(&args, bed_hash, &av, 1);
    } else {
        ret = amplicon_stats(&args, bed_hash, &argv[optind], argc-optind);
    }

    free(args.argv);
    destroy_bed_hash(bed_hash);

    return ret;
}
