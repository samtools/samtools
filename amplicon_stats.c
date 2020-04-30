/*  stats.c -- This is the former bamcheck integrated into samtools/htslib.

    Copyright (C) 2020 Genome Research Ltd.

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
 * - Cope with arbitrary sizes of reference and numbers of amplicons
 * - Cope with multiple references.  What do we do here?  Just request one?
 * - Permit regions rather than consuming whole file (maybe solves above).
 * - Detect double-length products (or primer pool hopping).
 * - eg 31_LEFT to 33_RIGHT, or even 31_LEFT to 32_RIGHT.
 *   Binning via TLEN field will suffice.
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

#include "samtools.h"
#include "sam_opts.h"
#include "bam_ampliconclip.h"

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(a) ((a)>=0?(a):-(a))
#endif

#define MAX_AMP_LEN 1000   // Maximum length of any single amplicon
#define MAX_PRIMER_PER_AMPLICON 4  // Max primers per LEFT/RIGHT
#define MAX_DEPTH 5        // Number of different depths permitted
#define MAX_TLEN 2000      // TLEN distribution size

typedef struct {
    sam_global_args ga;
    uint32_t flag_require;
    uint32_t flag_filter;
    bed_pair_list_t sites;
    int max_delta; // Used for matching read to amplicon primer loc
    int min_depth[MAX_DEPTH]; // Used for coverage; must be >= min_depth deep
    int use_sample_name;
    int max_amp;
    int64_t max_len;
    int tlen_adj; // Adjust tlen by this amount, due to clip but no fixmate
    FILE *out_fp;
} astats_args_t;

// We can have multiple primers for LEFT / RIGHT, so this
// permits detection by any compatible combination.
typedef struct {
    int64_t left[MAX_PRIMER_PER_AMPLICON];
    int nleft;
    int64_t right[MAX_PRIMER_PER_AMPLICON];
    int nright;
    int64_t max_left, min_right; // inner dimensions
    int64_t min_left, max_right; // outer dimensions
} amplicon_t;

// Map positions to amplicon numbers.
static int *pos2start = NULL;
static int *pos2end = NULL;
static int64_t pos_lookup_len = 0;

// Lookup table to go from position to amplicon based on
// read start / end.
//
// NB: Could do bed2amplicon and this code as bed2pos() in a single step.
static int initialise_amp_pos_lookup(astats_args_t *args,
                                     amplicon_t *amp, int namp,
                                     int64_t max_len) {
    int64_t i, j;

    if (!pos_lookup_len) {
        pos_lookup_len = max_len;
        if (!(pos2start = calloc(max_len+1, sizeof(*pos2start))))
            return -1;
        if (!(pos2end   = calloc(max_len+1, sizeof(*pos2start)))) {
            return -1;
        }
    } else if (pos_lookup_len != max_len) {
        fprintf(stderr, "[ampliconstats] error: "
                "input files with differing reference length");
        return -1;
    } else {
        return 0; // already done
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

static void free_amp_pos_lookup(void) {
    free(pos2start);
    free(pos2end);
}

// We're only interest in the internal part of the amplicon.
// Our bed file has LEFT start/end followed by RIGHT start/end,
// so collapse these to LEFT end / RIGHT start.
//
// Returns right most amplicon position on success,
//         < 0 on error
static int64_t bed2amplicon(astats_args_t *args,
                            bed_pair_list_t *sites,
                            amplicon_t *amp, int *namp) {
    int i, j;
    int64_t max_right = 0;

    *namp = 0;

    // Assume all primers for the same amplicon are adjacent in BED
    // with all + followed by all -.  Thus - to + signifies next primer set.
    int last_rev = 0;
    amp[0].max_left = 0;
    amp[0].min_right = INT64_MAX;
    amp[0].min_left = INT64_MAX;
    amp[0].max_right = 0;
    for (i = j = 0; i < sites->length; i++) {
        if (sites->bp[i].rev == 0 && last_rev) {
            j++;
            if (j >= args->max_amp) {
                fprintf(stderr, "[ampliconstats] error: "
                        "too many amplicons. Use -a option to raise this.");
                return -1;
            }
            amp[j].max_left = 0;
            amp[j].min_right = INT64_MAX;
            amp[j].min_left = INT64_MAX;
            amp[j].max_right = 0;
        }
        if (sites->bp[i].rev == 0) {
            if (amp[j].nleft >= MAX_PRIMER_PER_AMPLICON) {
                print_error_errno("ampliconstats",
                                  "too many primers per amplicon");
                return -1;
            }
            amp[j].left[amp[j].nleft++] = sites->bp[i].right;
            if (amp[j].max_left < sites->bp[i].right+1)
                amp[j].max_left = sites->bp[i].right+1;
            if (amp[j].min_left > sites->bp[i].right+1)
                amp[j].min_left = sites->bp[i].right+1;
        } else {
            if (amp[j].nright >= MAX_PRIMER_PER_AMPLICON) {
                print_error_errno("ampliconstats",
                                  "too many primers per amplicon");
                return -1;
            }
            amp[j].right[amp[j].nright++] = sites->bp[i].left;
            if (amp[j].min_right > sites->bp[i].left-1)
                amp[j].min_right = sites->bp[i].left-1;
            if (amp[j].max_right < sites->bp[i].left-1) {
                amp[j].max_right = sites->bp[i].left-1;
                if (amp[j].max_right - amp[j].min_left + 1 >= MAX_AMP_LEN) {
                    fprintf(stderr, "[ampliconstats] error: "
                            "amplicon longer than MAX_AMP_LEN define\n");
                    return -1;
                }
                if (max_right < amp[j].max_right)
                    max_right = amp[j].max_right;
            }
        }
        last_rev = sites->bp[i].rev;
    }
    *namp = ++j;

    if (j >= args->max_amp) {
        fprintf(stderr, "[ampliconstats] error: "
                "too many amplicons. Use -a option to raise this.");
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

typedef struct {
    int nseq;       // total sequence count
    int nfiltered;  // sequence filtered
    int nfailprimer;// count of sequences not matching the primer locations

    // Sizes of memory allocated below, to permit reset
    int max_amp;

    // Summary across all samples, sum(x) plus sum(x^2) for s.d. calc
    int64_t *nreads, *nreads2;          // [MAX_AMP]
    double  *nrperc, *nrperc2;          // [MAX_AMP]
    int64_t *nbases, *nbases2;          // [MAX_AMP]
    int64_t (*coverage)[MAX_AMP_LEN];   // [MAX_AMP][MAX_AMP_LEN]
    double  (*covered_perc)[MAX_DEPTH]; // [MAX_AMP][MAX_DEPTH]
    double  (*covered_perc2)[MAX_DEPTH];// [MAX_AMP][MAX_DEPTH];
    int     (*tlen_dist)[MAX_TLEN];     // [MAX_AMP][MAX_TLEN];

    // 0 is correct pair, 1 is incorrect pair, 2 is unidentified
    int     (*amp_dist)[3];             // [MAX_AMP][3];
    //int     amp_pair[MAX_AMP][MAX_AMP]; // dotplot style
} astats_t;

void stats_free(astats_t *st) {
    if (!st)
        return;

    free(st->nreads);
    free(st->nreads2);
    free(st->nrperc);
    free(st->nrperc2);
    free(st->nbases);
    free(st->nbases2);
    free(st->coverage);
    free(st->covered_perc);
    free(st->covered_perc2);
    free(st->tlen_dist);
    free(st->amp_dist);

    free(st);
}

astats_t *stats_alloc(int64_t max_len, int max_amp) {
    astats_t *st = calloc(1, sizeof(*st));
    if (!st)
        return NULL;

    st->max_amp = max_amp;

    if (!(st->nreads  = calloc(max_amp, sizeof(*st->nreads))))  goto err;
    if (!(st->nreads2 = calloc(max_amp, sizeof(*st->nreads2)))) goto err;
    if (!(st->nrperc  = calloc(max_amp, sizeof(*st->nrperc))))  goto err;
    if (!(st->nrperc2 = calloc(max_amp, sizeof(*st->nrperc2)))) goto err;
    if (!(st->nbases  = calloc(max_amp, sizeof(*st->nbases))))  goto err;
    if (!(st->nbases2 = calloc(max_amp, sizeof(*st->nbases2)))) goto err;

    if (!(st->coverage = calloc(max_amp, sizeof(*st->coverage)))) goto err;

    if (!(st->covered_perc  = calloc(max_amp, sizeof(*st->covered_perc))))
        goto err;
    if (!(st->covered_perc2 = calloc(max_amp, sizeof(*st->covered_perc2))))
        goto err;

    if (!(st->tlen_dist = calloc(max_amp, sizeof(*st->tlen_dist)))) goto err;
    if (!(st->amp_dist  = calloc(max_amp, sizeof(*st->amp_dist))))  goto err;

    return st;

 err:
    stats_free(st);
    return NULL;
}

void stats_reset(astats_t *st) {
    st->nseq = 0;
    st->nfiltered = 0;
    st->nfailprimer = 0;

    memset(st->nreads,  0, st->max_amp*sizeof(*st->nreads));
    memset(st->nreads2, 0, st->max_amp*sizeof(*st->nreads2));

    memset(st->nrperc,  0, st->max_amp*sizeof(*st->nrperc));
    memset(st->nrperc2, 0, st->max_amp*sizeof(*st->nrperc2));

    memset(st->nbases,  0, st->max_amp*sizeof(*st->nbases));
    memset(st->nbases2, 0, st->max_amp*sizeof(*st->nbases2));

    memset(st->coverage, 0, st->max_amp * sizeof(*st->coverage));
    memset(st->covered_perc,  0, st->max_amp * sizeof(*st->covered_perc));
    memset(st->covered_perc2, 0, st->max_amp * sizeof(*st->covered_perc2));

    memset(st->tlen_dist, 0, st->max_amp * sizeof(*st->tlen_dist));
    memset(st->amp_dist,  0, st->max_amp * sizeof(*st->amp_dist));
}

#define ISIZE_TO_DIST(i) (ABS(i)<MAX_TLEN?ABS(i):MAX_TLEN)

static int accumulate_stats(astats_t *stats,
                            astats_args_t *args,
                            amplicon_t *amp, int namp,
                            bam1_t *b) {
    stats->nseq++;
    if ((b->core.flag & args->flag_require) != args->flag_require ||
        (b->core.flag & args->flag_filter)  != 0) {
        stats->nfiltered++;
        return 0;
    }

    int64_t start = b->core.pos;
    int64_t end = bam_endpos(b);

    int anum = (b->core.flag & BAM_FREVERSE)
        ? (end-1 >= 0 && end-1 < args->max_len ? pos2end[end-1] : -1)
        : (start >= 0 && start < args->max_len ? pos2start[start] : -1);


    // ivar sometimes soft-clips 100% of the bases.
    // This is essentially unmapped
    if (end == start && (args->flag_filter & BAM_FUNMAP)) {
        stats->nfiltered++;
        return 0;
    }

    if (anum == -1) {
        stats->nfailprimer++;
        return 0;
    }

    stats->nreads[anum]++;
    stats->nbases[anum] += end-start; // NB: ref bases rather than read bases

    int64_t i;
    if (start < 0) start = 0;
    if (end > args->max_len) end = args->max_len;

    int64_t ostart = MAX(start, amp[anum].min_left-1);
    int64_t oend = MIN(end, amp[anum].max_right);
    int64_t offset = amp[anum].min_left-1;
    for (i = ostart; i < oend; i++)
        stats->coverage[anum][i-offset]++;

    // Basic template length
    stats->tlen_dist[anum][ISIZE_TO_DIST(b->core.isize)]++;

    // Template length in terms of amplicon number to amplicon number.
    // We expect left to right of same amplicon (len 0), but it may go
    // to next amplicon (len 1) or prev (len -1), etc.
    int64_t t_end = (b->core.flag & BAM_FREVERSE ? end : start)
        + b->core.isize;

    // If we've clipped the primers but not followed up with a fixmates
    // then our start+TLEN will take us to a location which is
    // length(LEFT_PRIMER) + length(RIGHT_PRIMER) too far away.
    //
    // The correct solution is to run samtools fixmate so TLEN is correct.
    // The hacky solution is to fudge the expected tlen by double the
    // average primer length (e.g. 50).
    t_end += b->core.isize > 0 ? -args->tlen_adj : +args->tlen_adj;
    int oth_anum = -1;
    if (t_end > 0 && t_end < args->max_len && b->core.isize != 0)
        oth_anum = (b->core.flag & BAM_FREVERSE)
            ? pos2start[t_end]
            : pos2end[t_end];

    // We don't want to count our pairs twice.
    // If both left/right are known, count it on left only.
    // If only one is known, we'll only get to this code once
    // so we can also count it.

    if (anum != -1 && oth_anum != -1) {
        if (start <= t_end)
            stats->amp_dist[anum][oth_anum == anum ? 0 : 1]++;
    } else {
        stats->amp_dist[anum][2]++;
    }

    return 0;
}

// Append file local stats to global stats
void append_stats(astats_t *lstats, astats_t *gstats, int namp) {
    gstats->nseq += lstats->nseq;
    gstats->nfiltered += lstats->nfiltered;
    gstats->nfailprimer += lstats->nfailprimer;

    int a;
    for (a = 0; a < namp; a++) {
        gstats->nreads[a]  += lstats->nreads[a];
        gstats->nreads2[a] += lstats->nreads[a] * lstats->nreads[a];

        // To get mean & sd for amplicon read percentage, we need
        // to do the divisions here as nseq differs for each sample.
        int nseq = lstats->nseq - lstats->nfiltered - lstats->nfailprimer;
        double nrperc = nseq ? 100.0 * lstats->nreads[a] / nseq : 0;
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

        for (d = 0; d < MAX_TLEN; d++)
            gstats->tlen_dist[a][d] += lstats->tlen_dist[a][d];

        for (d = 0; d < 3; d++)
            gstats->amp_dist[a][d] += lstats->amp_dist[a][d];
    }
}

void dump_stats(char type, char *name, astats_t *stats, astats_args_t *args,
                amplicon_t *amp, int namp, int nfile) {
    int i, j;
    FILE *ofp = args->out_fp;

    // summary stats for this sample (or for all samples)
    fprintf(ofp, "# Summary stats.\n");
    fprintf(ofp, "# Use 'grep ^%cSS | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cSS\t%s\traw total sequences:\t%d\n",
            type, name, stats->nseq);
    fprintf(ofp, "%cSS\t%s\tfiltered sequences:\t%d\n",
            type, name, stats->nfiltered);
    fprintf(ofp, "%cSS\t%s\tfailed primer match:\t%d\n",
            type, name, stats->nfailprimer);
    int nmatch = stats->nseq - stats->nfiltered - stats->nfailprimer;
    fprintf(ofp, "%cSS\t%s\tmatching sequences:\t%d\n",
            type, name, nmatch);

    int d = 0;
    do {
        // From first to last amplicon only, so not entire consensus.
        // If contig length is known, maybe we want to add the missing
        // count to < DEPTH figures?
        int64_t start = 0, covered = 0, total = 0;
        for (i = 0; i < namp; i++) {
            int64_t j, offset = amp[i].min_left-1;
            for (j = MAX(start, amp[i].max_left-1);
                 j < MAX(start, amp[i].min_right); j++) {
                if (stats->coverage[i][j-offset] >= args->min_depth[d])
                    covered++;
                total++;
            }
            start = MAX(start, amp[i].min_right);
        }
        fprintf(ofp, "%cSS\t%s\tconsensus depth count < %d and >= %d:\t%"
                PRId64"\t%"PRId64"\n", type, name,
                args->min_depth[d], args->min_depth[d],
                total-covered, covered);
    } while (++d < MAX_DEPTH && args->min_depth[d]);

    // Read count
    fprintf(ofp, "# Absolute matching read counts per amplicon.\n");
    fprintf(ofp, "# Use 'grep ^%cREADS | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cREADS\t%s", type, name);
    for (i = 0; i < namp; i++) {
        fprintf(ofp, "\t%"PRId64, stats->nreads[i]);
    }
    fprintf(ofp, "\n");

    if (type == 'C') {
        // For combined we can compute mean & standard deviation too
        fprintf(ofp, "CREADS\tMEAN");
        for (i = 0; i < namp; i++) {
            fprintf(ofp, "\t%.1f", stats->nreads[i] / (double)nfile);
        }
        fprintf(ofp, "\n");

        fprintf(ofp, "CREADS\tSTDDEV");
        for (i = 0; i < namp; i++) {
            double n1 = stats->nreads[i];
            fprintf(ofp, "\t%.1f", sqrt(stats->nreads2[i]/nfile - (n1/nfile)*(n1/nfile)));
        }
        fprintf(ofp, "\n");
    }

    fprintf(ofp, "# Read percentage of distribution between amplicons.\n");
    fprintf(ofp, "# Use 'grep ^%cRPERC | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cRPERC\t%s", type, name);
    for (i = 0; i < namp; i++) {
        if (type == 'C') {
            fprintf(ofp, "\t%.3f", (double)stats->nrperc[i] / nfile);
        } else {
            int nseq = stats->nseq - stats->nfiltered - stats->nfailprimer;
            fprintf(ofp, "\t%.3f", nseq ? 100.0 * stats->nreads[i] / nseq : 0);
        }
    }
    fprintf(ofp, "\n");

    if (type == 'C') {
        // For combined we compute mean and standard deviation too
        fprintf(ofp, "CRPERC\tMEAN");
        for (i = 0; i < namp; i++) {
            fprintf(ofp, "\t%.3f", stats->nrperc[i] / nfile);
        }
        fprintf(ofp, "\n");

        fprintf(ofp, "CRPERC\tSTDDEV");
        for (i = 0; i < namp; i++) {
            // variance = SUM(X^2) - ((SUM(X)^2) / N)
            double n1 = stats->nrperc[i];
            double v = stats->nrperc2[i]/nfile - (n1/nfile)*(n1/nfile);
            fprintf(ofp, "\t%.3f", sqrt(v));
        }
        fprintf(ofp, "\n");
    }

    // Base depth
    fprintf(ofp, "# Read depth per amplicon.\n");
    fprintf(ofp, "# Use 'grep ^%cREADS | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cDEPTH\t%s", type, name);
    for (i = 0; i < namp; i++) {
        int nseq = stats->nseq - stats->nfiltered - stats->nfailprimer;
        int64_t alen = amp[i].min_right - amp[i].max_left+1;
        fprintf(ofp, "\t%.1f", nseq ? stats->nbases[i] / (double)alen : 0);
    }
    fprintf(ofp, "\n");

    if (type == 'C') {
        // For combined we can compute mean & standard deviation too
        int nseq = stats->nseq - stats->nfiltered - stats->nfailprimer;
        fprintf(ofp, "CDEPTH\tMEAN");
        for (i = 0; i < namp; i++) {
            int64_t alen = amp[i].min_right - amp[i].max_left+1;
            fprintf(ofp, "\t%.1f", nseq ? stats->nbases[i] / (double)alen / nfile : 0);
        }
        fprintf(ofp, "\n");

        fprintf(ofp, "CDEPTH\tSTDDEV");
        for (i = 0; i < namp; i++) {
            double alen = amp[i].min_right - amp[i].max_left+1;
            double n1 = stats->nbases[i] / alen;
            double v = stats->nbases2[i] / (alen*alen) /nfile
                - (n1/nfile)*(n1/nfile);
            fprintf(ofp, "\t%.1f", sqrt(v));
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
            for (i = 0; i < namp; i++) {
                int covered = 0;
                int64_t j, offset = amp[i].min_left-1;
                for (j = amp[i].max_left-1; j < amp[i].min_right; j++)
                    if (stats->coverage[i][j-offset] >= args->min_depth[d])
                        covered++;
                int64_t alen = amp[i].min_right - amp[i].max_left+1;
                stats->covered_perc[i][d] = 100.0 * covered / alen;
                fprintf(ofp, "\t%.2f", 100.0 * covered / alen);
            }
            fprintf(ofp, "\n");
        } while (++d < MAX_DEPTH && args->min_depth[d]);
    } else if (type == 'C') {
        // For combined we can compute mean & standard deviation too
        int d = 0;
        do {
            fprintf(ofp, "CPCOV-%d\tMEAN", args->min_depth[d]);
            for (i = 0; i < namp; i++) {
                fprintf(ofp, "\t%.1f", stats->covered_perc[i][d] / nfile);
            }
            fprintf(ofp, "\n");

            fprintf(ofp, "CPCOV-%d\tSTDDEV", args->min_depth[d]);
            for (i = 0; i < namp; i++) {
                double n1 = stats->covered_perc[i][d] / nfile;
                double v = stats->covered_perc2[i][d] / nfile - n1*n1;
                fprintf(ofp, "\t%.1f", sqrt(v));
            }
            fprintf(ofp, "\n");
        } while (++d < MAX_DEPTH && args->min_depth[d]);
    }


    // TLEN distribution
    fprintf(ofp, "# Distribution of aligned template length.\n");
    fprintf(ofp, "# Use 'grep ^%cTLEN | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cTLEN\t%s\t0", type, name); // all merged
    int tlen_dist[MAX_TLEN] = {0};
    for (i = 0; i < namp; i++) // accumulate for all amps
        for (j = 0; j < MAX_TLEN; j++)
            tlen_dist[j] += stats->tlen_dist[i][j];
    for (j = 0; j < MAX_TLEN; j++)
        fprintf(ofp, "\t%d", tlen_dist[j]);
    fprintf(ofp, "\n");

    for (i = 0; i < namp; i++) {
        fprintf(ofp, "%cTLEN\t%s\t%d", type, name, i+1); // per amplicon
        for (j = 0; j < MAX_TLEN; j++)
            fprintf(ofp, "\t%d", stats->tlen_dist[i][j]);
        fprintf(ofp, "\n");
    }


    // AMP length distribution.
    // 0 = both ends in this amplicon
    // 1 = ends in different amplicons
    // 2 = other end matching an unknown amplicon site
    //     (see tlen_dist for further analysis of where)
    fprintf(ofp, "# Classification of amplicon status.  Columns are\n");
    fprintf(ofp, "# number with both primers from this amplicon, number with\n");
    fprintf(ofp, "# primers from different amplicon, and number with a position\n");
    fprintf(ofp, "# not matching any valid amplicon primer site\n");
    fprintf(ofp, "# Use 'grep ^%cAMP | cut -f 2-' to extract this part.\n", type);
    fprintf(ofp, "%cAMP\t%s\t0", type, name); // all merged
    int amp_dist[3] = {0};
    for (i = 0; i < namp; i++) { // accumulate for all amps
        amp_dist[0] += stats->amp_dist[i][0];
        amp_dist[1] += stats->amp_dist[i][1];
        amp_dist[2] += stats->amp_dist[i][2];
    }
    fprintf(ofp, "\t%d\t%d\t%d\n", amp_dist[0], amp_dist[1], amp_dist[2]);

    for (i = 0; i < namp; i++) {
        fprintf(ofp, "%cAMP\t%s\t%d", type, name, i+1); // per amplicon
        fprintf(ofp, "\t%d\t%d\t%d\n", stats->amp_dist[i][0],
                stats->amp_dist[i][1], stats->amp_dist[i][2]);
    }

//    for (i = 0; i < namp; i++) {
//        printf("%cAMP\t%s\t%d", type, name, i+1); // per amplicon
//        for (j = 0; j < namp; j++)
//            printf("\t%d", stats->amp_pair[i][j]);
//        printf("\n");
//    }

}

char const *get_sample_name(sam_hdr_t *header, char *RG) {
    kstring_t ks = {0};
    sam_hdr_find_tag_id(header, "RG", RG?"ID":NULL, RG, "SM", &ks);
    return ks.s;
}

int64_t get_ref_len(sam_hdr_t *header, char *SQ) {
    int tid = SQ ? sam_hdr_name2tid(header, SQ) : 0;
    return tid >= 0 ? sam_hdr_tid2len(header, tid) : -1;
}

static int amplicon_stats(astats_args_t *args, char **filev, int filec) {
    int i, namp;
    samFile *fp = NULL;
    sam_hdr_t *header = NULL;
    bam1_t *b = bam_init1();
    FILE *ofp = args->out_fp;
    char sname_[8192], *sname = NULL;

    // Global stats across all samples
    astats_t *gstats = NULL, *lstats = NULL;

    amplicon_t *amp = calloc(args->sites.length, sizeof(*amp));
    if (!amp)
        return -1;

    if (bed2amplicon(args, &args->sites, amp, &namp) < 0)
        goto err;

    fprintf(ofp, "# Summary statistics, used for scaling the plots.\n");
    fprintf(ofp, "SS\tNumber of amplicons:\t%d\n", namp);
    fprintf(ofp, "SS\tNumber of files:\t%d\n", filec);
    fprintf(ofp, "SS\tEnd of summary\n");

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

        if (args->use_sample_name)
            sname = (char *)get_sample_name(header, NULL);

        // FIXME: permit other references to be specified.
        if ((args->max_len = get_ref_len(header, NULL)) < 0)
            goto err;
        if (initialise_amp_pos_lookup(args, amp, namp, args->max_len) < 0)
            goto err;

        if (!gstats) gstats = stats_alloc(args->max_len, args->max_amp);
        if (!lstats) lstats = stats_alloc(args->max_len, args->max_amp);
        if (!lstats || !gstats)
            goto err;


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
        stats_reset(lstats);

        int r;
        while ((r = sam_read1(fp, header, b)) >= 0) {
            if (accumulate_stats(lstats, args, amp, namp, b) < 0)
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

        dump_stats('F', sname, lstats, args, amp, namp, filec);

        append_stats(lstats, gstats, namp);

        if (sname && sname != sname_)
            free(sname);
        sname = NULL;
    }

    dump_stats('C', "COMBINED", gstats, args, amp, namp, filec);

    stats_free(lstats);
    stats_free(gstats);
    bam_destroy1(b);
    free(amp);
    free_amp_pos_lookup();

    return 0;

 err:
    bam_destroy1(b);
    if (header)
        sam_hdr_destroy(header);
    if (fp)
        sam_close(fp);
    free(amp);
    free_amp_pos_lookup();
    if (sname && sname != sname_)
        free(sname);

    return -1;
}

static int usage(astats_args_t *args, FILE *fp, int exit_status) {
    fprintf(fp,
"\n"
"Usage: samtools ampliconstats [options] primers.bed *.bam > astats.txt\n"
"\n"
"Options:\n");
    fprintf(fp, "  -f, --required-flag STR|INT\n"
            "               Only include reads with all  of the FLAGs in present [0x%X]\n",args->flag_require);
    fprintf(fp, "  -F, --filter-flag STR|INT\n"
            "               Only include reads with none of the FLAGs in present [0x%X]\n",args->flag_filter & 0xffff);
    fprintf(fp, "  -a, --max-amplicons INT\n"
            "               Change the maximum number of amplicons permitted\n");
    fprintf(fp, "  -d, --min-depth INT[,INT]...\n"
            "               Minimum base depth(s) to consider position covered [%d]\n", args->min_depth[0]);
    fprintf(fp, "  -m, --pos-margin INT\n"
            "               Margin of error for matching primer positions [%d].\n", args->max_delta);
    fprintf(fp, "  -o, FILE\n"
            "               Specify output file [stdout if unset]\n");
    fprintf(fp, "  -s, --use-sample-name\n"
            "               Use the sample name from the first @RG header line\n");
    fprintf(fp, "  -t, --tlen-adjust INT\n"
            "               Add/subtract from TLEN; use when clipping but no fixmate step\n");
    sam_global_opt_help(fp, "I.--.@");

    return exit_status;
}

int main_ampliconstats(int argc, char **argv) {
    astats_args_t args = {
        .ga = SAM_GLOBAL_ARGS_INIT,
        .flag_require = 0,
        .flag_filter = 0x10F04,
        .sites = {NULL, 0, 0},
        .max_delta = 30, // large enough to cope with alt primers
        .min_depth = {1},
        .use_sample_name = 0,
        .max_amp = 1000,
        .tlen_adj = 0,
        .out_fp = stdout,
    }, oargs = args;

    static const struct option loptions[] =
    {
        SAM_OPT_GLOBAL_OPTIONS('I', 0, '-', '-', 0, '@'),
        {"help", no_argument, NULL, 'h'},
        {"flag-require", required_argument, NULL, 'f'},
        {"flag-filter", required_argument, NULL, 'F'},
        {"min-depth", required_argument, NULL, 'd'},
        {"pos-margin", required_argument, NULL, 'm'},
        {"use-sample-name", no_argument, NULL, 's'},
        {"max-amplicons", required_argument, NULL, 'a'},
        {"tlen-adjust", required_argument, NULL, 't'},
        {NULL, 0, NULL, 0}
    };
    int opt;

    while ( (opt=getopt_long(argc,argv,"?hf:F:@:p:m:d:sa:t:o:",loptions,NULL))>0 ) {
        switch (opt) {
        case 'f': args.flag_require = bam_str2flag(optarg); break;
        case 'F':
            if (args.flag_filter & 0x10000)
                args.flag_filter = 0; // strip default on first -F usage
            args.flag_filter |= bam_str2flag(optarg); break;

        case 'm': args.max_delta = atoi(optarg); break; // margin
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

        case 't': args.tlen_adj = atoi(optarg);break;

        case 's': args.use_sample_name = 1;break;

        case 'o':
            if (!(args.out_fp = fopen(optarg, "w"))) {
                perror(optarg);
                return 1;
            }
            break;

        case '?': return usage(&oargs, stderr, EXIT_FAILURE);
        case 'h': return usage(&oargs, stdout, EXIT_SUCCESS);

        default:
            if (parse_sam_global_opt(opt, optarg, loptions, &args.ga) != 0)
                usage(&oargs,stderr, EXIT_FAILURE);
            break;
        }
    }

    if (argc <= optind+1 && isatty(STDIN_FILENO))
        return usage(&oargs, stdout, EXIT_SUCCESS);

    int64_t longest;
    if (load_bed_file_pairs(argv[optind++], 1, 0, &args.sites, &longest)) {
        print_error_errno("ampliconstats",
                          "Could not read file \"%s\"", optarg);
        return 1;
    }

    int ret;
    if (argc == optind) {
        char *av = "-";
        ret = amplicon_stats(&args, &av, 1);
    } else {
        ret = amplicon_stats(&args, &argv[optind], argc-optind);
    }

    free(args.sites.bp);

    return ret;
}
