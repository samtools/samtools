/*  bam2depth.c -- depth subcommand.

    Copyright (C) 2011, 2012 Broad Institute.
    Copyright (C) 2012-2016, 2018, 2019-2022 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk> (to 2020)
    Author: James Bonfield <jkb@sanger.ac.uk> (2021 rewrite)


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

/* This program demonstrates how to generate pileup from multiple BAMs
 * simultaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -lhts -lz
 */

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "samtools.h"
#include "bedidx.h"
#include "sam_opts.h"
#include "htslib/khash.h"

// From bam_plcmd.c
int read_file_list(const char *file_list, int *n, char **argv[]);

// We accumulate to hist[pos & (size-1)].  This is a ring-buffer.
// We track where we last got to in output and what the biggest value
// we've written to so far (in absolute unmasked coordinates) in
// "last_output" and "end_pos" respectively.
// For each new record we just flush anything we haven't written yet
// already, between "last_output" and this read's start position, and
// initialise any newly seen positions between "end_pos" and this read's
// end position.
typedef struct {
    size_t size;
    int **hist;         // hist[nfiles][size]
    hts_pos_t *end_pos; // end_pos[nfiles]
    hts_pos_t last_output;
    int last_ref;
    int nfiles;
    const char *ref;
    kstring_t ks;
    hts_pos_t beg, end; // limit to region
    int tid;
} depth_hist;

typedef struct {
    int header;
    int flag;
    int incl_flag;
    int require_flag;
    int min_qual;
    int min_mqual;
    int min_len;
    int skip_del;
    int all_pos;
    int remove_overlaps;
    FILE *out;
    char *reg;
    void *bed;
} depth_opt;

static void zero_region(depth_opt *opt, depth_hist *dh,
                        const char *name, hts_pos_t start, hts_pos_t end) {
    hts_pos_t i;
    kstring_t *ks = &dh->ks;

    kputs(name, ks_clear(ks));
    kputc('\t', ks);
    size_t cur_l = ks->l;
    if (dh->beg >= 0 && start < dh->beg)
        start = dh->beg;
    if (dh->end >= 0 && end > dh->end)
        end = dh->end;

    for (i = start; i < end; i++) {
        // Could be optimised, but needs better API to skip to next
        // bed region.
        if (opt->bed && bed_overlap(opt->bed, name, i, i+1) == 0)
            continue;

        ks->l = cur_l;
        kputll(i+1,  ks);
        int n;
        for (n = 0; n < dh->nfiles; n++) {
            kputc_('\t', ks);
            kputc_('0',  ks);
        }
        kputc('\n',  ks);
        fputs(ks->s, opt->out);
    }
    ks->l = cur_l;
}

// A variation of bam_cigar2qlen which doesn't count soft-clips in to the
// equation.  Basically it's the number of bases in query that are aligned
// in some way to the reference (including insertions, which are considered
// to be aligned by dint of being anchored either side).
hts_pos_t qlen_used(bam1_t *b) {
    int n_cigar = b->core.n_cigar;
    const uint32_t *cigar = bam_get_cigar(b);

    hts_pos_t l;

    if (b->core.l_qseq) {
        // Known SEQ permits of short cut of l_qseq minus CSOFT_CLIPs.
        // Full scan not needed, which helps on excessively long CIGARs.
        l = b->core.l_qseq;
        int kl, kr;
        for (kl = 0; kl < n_cigar; kl++)
            if (bam_cigar_op(cigar[kl]) == BAM_CSOFT_CLIP)
                l -= bam_cigar_oplen(cigar[kl]);
            else
                break;

        for (kr = n_cigar-1; kr > kl; kr--)
            if (bam_cigar_op(cigar[kr]) == BAM_CSOFT_CLIP)
                l -= bam_cigar_oplen(cigar[kr]);
            else
                break;
    } else {
        // Unknown SEQ ("*") needs a full scan through the CIGAR string.
        static int query[16] = {
          //M I D N  S H P =  X B ? ?  ? ? ? ?
            1,1,0,0, 0,0,0,1, 1,0,0,0, 0,0,0,0
        };
        int k;
        for (k = l = 0; k < n_cigar; k++)
            if (query[bam_cigar_op(cigar[k])])
                l += bam_cigar_oplen(cigar[k]);
    }
    return l;

}

// Adds the depth for a single read to a depth_hist struct.
// For just one file, this is easy.  We just have a circular buffer
// where we increment values for bits that overlap existing data
// and initialise values for coordinates which we're seeing for the first
// time.  This is tracked by "end_pos" to know where we've got to.
//
// As the input is sorted, we can flush output from "last_output" to
// b->core.pos.
//
// With multiple files, we must feed data in sorted order as if all files
// are merged, but track depth per file.  This also means "end_pos" is per
// file too, but "last_output" is global as it corresponds to rows printed.
static int add_depth(depth_opt *opt, depth_hist *dh, sam_hdr_t *h, bam1_t *b,
                     int overlap_clip, int file) {
    hts_pos_t i;
    size_t hmask = dh->size-1;
    int n;

    if (!b || b->core.tid != dh->last_ref) {
        // New ref
        if (dh->last_ref >= 0) {
            // do end
            size_t cur_l = dh->ks.l;
            int nf = dh->nfiles;
            i = dh->last_output;
            for (i = dh->last_output; nf; i++) {
                nf = 0;
                for (n = 0; n < dh->nfiles; n++) {
                    if (i < dh->end_pos[n])
                        nf++;
                }
                if (!nf)
                    break;

                if (opt->bed && bed_overlap(opt->bed, dh->ref, i, i+1) == 0)
                    continue;

                dh->ks.l = cur_l;
                kputll(i+1, &dh->ks);
                for (n = 0; n < dh->nfiles; n++) {
                    kputc_('\t', &dh->ks);
                    int d = i < dh->end_pos[n]
                        ? dh->hist[n][i & hmask]
                        : 0;
                    kputuw(d, &dh->ks);
                }
                kputc('\n', &dh->ks);
                fputs(dh->ks.s, opt->out);
            }
            if (opt->all_pos) {
                // End of last ref
                zero_region(opt, dh,
                            sam_hdr_tid2name(h, dh->last_ref),
                            i, sam_hdr_tid2len(h, dh->last_ref));
            }
            dh->ks.l = cur_l;
        }

        if (opt->all_pos > 1 && !opt->reg) {
            // Any previous unused refs
            int lr = dh->last_ref < 0 ? 0 : dh->last_ref+1;
            int rr = b ? b->core.tid : sam_hdr_nref(h), r;
            for (r = lr; r < rr; r++)
                zero_region(opt, dh,
                            sam_hdr_tid2name(h, r),
                            0, sam_hdr_tid2len(h, r));
        }

        if (!b) {
            // we're just flushing to end of file
            if (opt->all_pos && opt->reg && dh->last_ref < 0)
                // -a or -aa without a single read being output yet
                zero_region(opt, dh, sam_hdr_tid2name(h, dh->tid), dh->beg,
                            MIN(dh->end, sam_hdr_tid2len(h, dh->tid)));

            return 0;
        }

        for (n = 0; dh->end_pos && n < dh->nfiles; n++)
            dh->end_pos[n] = 0;
        dh->last_output = dh->beg >= 0
            ? MAX(b->core.pos, dh->beg)
            : b->core.pos;
        dh->last_ref = b->core.tid;
        dh->ref = sam_hdr_tid2name(h, b->core.tid);
        kputs(dh->ref, ks_clear(&dh->ks));
        kputc('\t', &dh->ks);

        if (opt->all_pos)
            // Start of ref
            zero_region(opt, dh, dh->ref, 0, b->core.pos);
    } else {
        if (dh->last_output < b->core.pos) {
            // Flush any depth outputs up to start of new read
            size_t cur_l = dh->ks.l;
            int nf = dh->nfiles;
            for (i = dh->last_output; i < b->core.pos; i++) {
                nf = 0;
                for (n = 0; n < dh->nfiles; n++) {
                    if (i < dh->end_pos[n])
                        nf++;
                }
                if (!nf)
                    break;

                if (opt->bed && bed_overlap(opt->bed, dh->ref, i, i+1) == 0)
                    continue;

                dh->ks.l = cur_l;
                kputll(i+1, &dh->ks);
                for (n = 0; n < dh->nfiles; n++) {
                    kputc_('\t', &dh->ks);
                    int d = i < dh->end_pos[n]
                        ? dh->hist[n][i & hmask]
                        : 0;
                    kputuw(d, &dh->ks);
                }
                kputc('\n', &dh->ks);
                fputs(dh->ks.s, opt->out);
            }
            if (opt->all_pos && i < b->core.pos)
                // Hole in middle of ref
                zero_region(opt, dh, dh->ref, i, b->core.pos);

            dh->ks.l = cur_l;
            dh->last_output = b->core.pos;
        }
    }

    hts_pos_t end_pos = bam_endpos(b); // 0 based, 1 past end.
    //printf("%d %d\n", (int)b->core.pos+1, (int)end_pos);

    if (b->core.tid < dh->last_ref ||
        (dh->last_ref == b->core.tid && end_pos < dh->last_output)) {
        print_error_errno("depth", "Data is not position sorted");
        return -1;
    }

    // If needed, grow the circular buffer.
    if (end_pos+1 - b->core.pos >= dh->size) {
        size_t old_size = dh->size;
        size_t old_hmask = hmask;
        while (end_pos+1 - b->core.pos >= dh->size)
            dh->size = dh->size ? 2*dh->size : 2048;
        hmask = dh->size-1;
        if (!dh->hist) {
            dh->hist = calloc(dh->nfiles, sizeof(*dh->hist));
            dh->end_pos = calloc(dh->nfiles, sizeof(*dh->end_pos));
            if (!dh->hist || !dh->end_pos)
                return -1;
        }
        for (n = 0; n < dh->nfiles; n++) {
            int *hist = calloc(dh->size, sizeof(*dh->hist[n]));
            if (!hist)
                return -1;

            // Simple approach for now; copy over old histogram verbatim.
            for (i = dh->last_output; i < dh->last_output + old_size; i++)
                hist[i & hmask] = dh->hist[n][i & old_hmask];
            free(dh->hist[n]);
            dh->hist[n] = hist;
        }
    }

    // Accumulate depth, based on CIGAR
    uint32_t *cig = bam_get_cigar(b);
    int ncig = b->core.n_cigar, j, k, spos = 0;

    // Zero new (previously unseen) coordinates so increment works later.
    hts_pos_t end = MAX(dh->end_pos[file], b->core.pos);
    if (end_pos > end && (end & hmask) < (end_pos & hmask)) {
        memset(&dh->hist[file][end & hmask], 0,
               sizeof(**dh->hist) * (end_pos - end));
    } else {
        for (i = end; i < end_pos; i++)
            dh->hist[file][i & hmask] = 0;
    }

    i = b->core.pos;
    uint8_t *qual = bam_get_qual(b);
    int min_qual = opt->min_qual;
    for (j = 0; j < ncig; j++) {
        int op    = bam_cigar_op(cig[j]);
        int oplen = bam_cigar_oplen(cig[j]);

        switch (op) {
        case BAM_CDEL:
        case BAM_CREF_SKIP:
            if (op != BAM_CDEL || opt->skip_del) {
                // don't increment reference location
                if (i + oplen >= dh->end_pos[file]) {
                    for (k = 0; k < oplen; k++, i++) {
                        if (i >= dh->end_pos[file])
                            // redundant due to zero new elements above?
                            dh->hist[file][i & hmask] = 0;
                    }
                } else {
                    i += oplen;
                }
            } else { // op == BAM_CDEL and we count them (-J option),
                // We don't incr spos here, but we still use qual.
                // This doesn't make much sense, but it's for compatibility
                // with the old code.  Arguably DEL shouldn't have a min
                // qual and should always pass (as we've explicitly asked to
                // include them).
                int *hist = dh->hist[file];
                k = 0;
                if (overlap_clip) {
                    if (i+oplen < overlap_clip) {
                        i += oplen;
                        break;
                    } else if (i < overlap_clip) {
                        k = overlap_clip - i;
                        i = overlap_clip;
                    }
                }

                // Question: should we even check quality values for DEL?
                // We've explicitly asked to include them, and the quality
                // is wrong anyway (it's the neighbouring base).  We do this
                // for now for compatibility with the old depth command.

                if (spos < b->core.l_qseq)
                    for (; k < oplen; k++, i++)
                        hist[i & hmask]+=qual[spos]>=min_qual;
                else
                    for (; k < oplen; k++, i++)
                        hist[i & hmask]++;
            }
            break;

        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            if ((i & hmask) < ((i+oplen) & hmask)) {
                // Optimisation when not wrapping around

                // Unrolling doesn't help clang, but helps gcc,
                // especially when not using -O3.
                int *hist = &dh->hist[file][i & hmask];
                if (min_qual || overlap_clip) {
                    k = 0;
                    if (overlap_clip) {
                        if (i+oplen < overlap_clip) {
                            i += oplen;
                            spos += oplen;
                            break;
                        } else if (i < overlap_clip) {
                            oplen -= overlap_clip - i;
                            spos += overlap_clip - i;
                            hist += overlap_clip - i;
                            i = overlap_clip;
                        }
                    }

                    // approx 50% of this func cpu time in this loop
                    for (; k < (oplen & ~7); k+=8) {
                        hist[k+0]+=qual[spos+0]>=min_qual;
                        hist[k+1]+=qual[spos+1]>=min_qual;
                        hist[k+2]+=qual[spos+2]>=min_qual;
                        hist[k+3]+=qual[spos+3]>=min_qual;
                        hist[k+4]+=qual[spos+4]>=min_qual;
                        hist[k+5]+=qual[spos+5]>=min_qual;
                        hist[k+6]+=qual[spos+6]>=min_qual;
                        hist[k+7]+=qual[spos+7]>=min_qual;
                        spos += 8;
                    }
                } else {
                    // easier to vectorize when no min_qual
                    for (k = 0; k < (oplen & ~7); k+=8) {
                        hist[k+0]++;
                        hist[k+1]++;
                        hist[k+2]++;
                        hist[k+3]++;
                        hist[k+4]++;
                        hist[k+5]++;
                        hist[k+6]++;
                        hist[k+7]++;
                    }
                    spos += k;
                }
                for (; k < oplen && spos < b->core.l_qseq; k++, spos++)
                    hist[k]+=qual[spos]>=min_qual;
                for (; k < oplen; k++, spos++)
                    hist[k]++;
                i += oplen;
            } else {
                // Simple to understand case, but slower.
                // We use this only for reads with wrap-around.
                int *hist = dh->hist[file];
                k = 0;
                if (overlap_clip) {
                    if (i+oplen < overlap_clip) {
                        i += oplen;
                        break;
                    } else if (i < overlap_clip) {
                        oplen -= overlap_clip - i;
                        spos += overlap_clip - i;
                        i = overlap_clip;
                    }
                }
                for (; k < oplen && spos < b->core.l_qseq; k++, i++, spos++)
                    hist[i & hmask]+=qual[spos]>=min_qual;
                for (; k < oplen; k++, i++, spos++)
                    hist[i & hmask]++;
            }
            break;

        case BAM_CINS:
        case BAM_CSOFT_CLIP:
            spos += oplen;
            break;

        case BAM_CPAD:
        case BAM_CHARD_CLIP:
            // ignore
            break;

        default:
            print_error("depth", "Unsupported cigar op '%d'", op);
            return -1;
        }
    }

    if (dh->end >= 0 && end_pos > dh->end)
        end_pos = dh->end;
    if (dh->end_pos[file] < end_pos)
        dh->end_pos[file] = end_pos;

    return 0;
}

// Hash on name -> alignment end pos. This permits a naive overlap removal.
// Note it cannot analyse the overlapping sequence and qualities, so the
// interaction of basecalls/qualities and the -Q parameter cannot be
// applied here (unlike the full mpileup algorithm).
KHASH_MAP_INIT_STR(olap_hash, hts_pos_t)
typedef khash_t(olap_hash) olap_hash_t;

static int fastdepth_core(depth_opt *opt, uint32_t nfiles, char **fn,
                          samFile **fp, hts_itr_t **itr, sam_hdr_t **h) {
    int ret = -1, err = 1, i;
    olap_hash_t **overlaps = NULL;
    depth_hist dh = {0};

    // An array of bam structs, one per input file, to hold the next entry
    bam1_t **b = calloc(nfiles, sizeof(*b));
    int *finished = calloc(nfiles, sizeof(*finished)), to_go = nfiles;
    if (!b || !finished)
        goto err;

    for (i = 0; i < nfiles; i++)
        if (!(b[i] = bam_init1()))
            goto err;

    // Do we need one overlap hash per file? Or shared?
    if (opt->remove_overlaps) {
        if (!(overlaps = calloc(nfiles, sizeof(*overlaps))))
            return -1;
        for (i = 0; i < nfiles; i++) {
            if (!(overlaps[i] = kh_init(olap_hash)))
                return -1;
        }
    }

    // Create the initial histogram
    dh.nfiles = nfiles;
    dh.size = 0;
    dh.hist = NULL;
    dh.last_ref = -99;
    dh.end_pos = NULL;
    dh.last_output = itr && itr[0] ? itr[0]->beg : 0;
    ks_initialize(&dh.ks);

    // Clip results to region if specified
    dh.beg = -1;
    dh.end = -1;
    dh.tid = 0;
    if (itr && itr[0]) {
        dh.tid = itr[0]->tid;
        dh.beg = itr[0]->beg;
        dh.end = itr[0]->end;
    }

    if (opt->header) {
        fprintf(opt->out, "#CHROM\tPOS");
        for (i = 0; i < nfiles; i++)
            fprintf(opt->out, "\t%s", fn[i]);
        fputc('\n', opt->out);
    }

    // Populate first record per file
    for (i = 0; i < nfiles; i++) {
        for(;;) {
            ret = itr && itr[i]
                ? sam_itr_next(fp[i], itr[i], b[i])
                : sam_read1(fp[i], h[i], b[i]);
            if (ret < -1)
                goto err;
            if (ret == -1) {
                to_go--;
                finished[i] = 1;
                break;
            }

            if (b[i]->core.tid < 0)
                continue;
            if (b[i]->core.flag & opt->flag)
                continue; // must have none of the flags set
            if (opt->incl_flag && (b[i]->core.flag & opt->incl_flag) == 0)
                continue; // must have at least one flag set
            if ((b[i]->core.flag & opt->require_flag) != opt->require_flag)
                continue; // must have all lags set
            if (b[i]->core.qual < opt->min_mqual)
                continue;

            // Original samtools depth used the total sequence (l_qseq)
            // including soft-clips.  This doesn't feel like a useful metric
            // to be filtering on.  We now only count sequence bases that
            // form the used part of the alignment.
            if (opt->min_len) {
                if (qlen_used(b[i]) < opt->min_len)
                    continue;
            }

            break;
        }
    }

    // Loop through input files, merging in order so we're
    // always adding the next record in sequence
    while (to_go) {
        // Find next record in file list
        int best_tid = INT_MAX, best_file = 0;
        hts_pos_t best_pos = HTS_POS_MAX;

        for (i = 0; i < nfiles; i++) {
            if (finished[i])
                continue;
            if (best_tid > b[i]->core.tid) {
                best_tid = b[i]->core.tid;
                best_pos = b[i]->core.pos;
                best_file = i;
            } else if (best_tid == b[i]->core.tid &&
                       best_pos > b[i]->core.pos) {
                best_pos = b[i]->core.pos;
                best_file = i;
            }
        }
        i = best_file;

        hts_pos_t clip = 0;
        if (overlaps && (b[i]->core.flag & BAM_FPAIRED) &&
            !(b[i]->core.flag & BAM_FMUNMAP)) {
            khiter_t k = kh_get(olap_hash, overlaps[i], bam_get_qname(b[i]));
            if (k == kh_end(overlaps[i])) {
                // not seen before
                hts_pos_t endpos = bam_endpos(b[i]);

                // Don't add if mate location is known and can't overlap.
                if (b[i]->core.mpos == -1 ||
                    (b[i]->core.tid == b[i]->core.mtid &&
                     b[i]->core.mpos <= endpos)) {
                    k = kh_put(olap_hash, overlaps[i], bam_get_qname(b[i]),
                               &ret);
                    if (ret < 0)
                        return -1;
                    kh_key(overlaps[i], k) = strdup(bam_get_qname(b[i]));
                    kh_value(overlaps[i], k) = endpos;
                }
            } else {
                // seen before
                clip = kh_value(overlaps[i], k);
                free((char *)kh_key(overlaps[i], k));
                kh_del(olap_hash, overlaps[i], k);
            }
        }

        // Add the next merged BAM record to the depth plot
        if ((ret = add_depth(opt, &dh, h[i], b[i], clip, i)) < 0) {
            ret = -1;
            goto err;
        }

        // Populate next record from this file
        for(;!finished[i];) {
            ret = itr && itr[i]
                ? sam_itr_next(fp[i], itr[i], b[i])
                : sam_read1(fp[i], h[i], b[i]);
            if (ret < -1) {
                ret = -1;
                goto err;
            }
            if (ret == -1) {
                to_go--;
                finished[i] = 1;
                break;
            }

            if (b[i]->core.tid < 0)
                continue;
            if (b[i]->core.flag & opt->flag)
                continue; // must have none of the flags set
            if (opt->incl_flag && (b[i]->core.flag & opt->incl_flag) == 0)
                continue; // must have at least one flag set
            if ((b[i]->core.flag & opt->require_flag) != opt->require_flag)
                continue; // must have all lags set
            if (b[i]->core.qual < opt->min_mqual)
                continue;

            if (opt->min_len) {
                if (qlen_used(b[i]) < opt->min_len)
                    continue;
            }

            break;
        }
    }

    // Tidy up end.
    ret = add_depth(opt, &dh, h[0], NULL, 0, 0);
    err = 0;

 err:
    if (ret == 0 && err)
        ret = -1;

    for (i = 0; i < nfiles; i++) {
        if (b[i])
            bam_destroy1(b[i]);
        if (dh.hist && dh.hist[i])
            free(dh.hist[i]);
    }
    free(b);
    free(finished);
    ks_free(&dh.ks);
    free(dh.hist);
    free(dh.end_pos);
    if (overlaps) {
        khiter_t k;
        for (i = 0; i < nfiles; i++) {
            if (!overlaps[i])
                continue;
            for (k = kh_begin(overlaps[i]); k < kh_end(overlaps[i]); k++)
                if (kh_exist(overlaps[i], k))
                    free((char *)kh_key(overlaps[i], k));
            kh_destroy(olap_hash, overlaps[i]);
        }
        free(overlaps);
    }

    return ret;
}

static void usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools depth [options] in.bam [in.bam ...]\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "  -a           Output all positions (including zero depth)\n");
    fprintf(fp, "  -a -a, -aa   Output absolutely all positions, including unused ref seqs\n");
    fprintf(fp, "  -r REG       Specify a region in chr or chr:from-to syntax\n");
    fprintf(fp, "  -b FILE      Use bed FILE for list of regions\n");
    fprintf(fp, "  -f FILE      Specify list of input BAM/SAM/CRAM filenames\n");
    fprintf(fp, "  -X           Use custom index files (in -X *.bam *.bam.bai order)\n");
    fprintf(fp, "  -g INT       Remove specified flags from default filter-out flag list\n");
    fprintf(fp, "  -G, --excl-flags FLAGS\n");
    fprintf(fp, "               Add specified flags to the  default filter-out flag list\n");
    fprintf(fp, "               [UNMAP,SECONDARY,QCFAIL,DUP]\n");
    fprintf(fp, "      --incl-flags FLAGS\n");
    fprintf(fp, "               Only include records with at least one the FLAGs present [0]\n");
    fprintf(fp, "      --require-flags FLAGS\n");
    fprintf(fp, "               Only include records with all of the FLAGs present [0]\n");
    fprintf(fp, "  -H           Print a file header line\n");
    fprintf(fp, "  -l INT       Minimum read length [0]\n");
    fprintf(fp, "  -o FILE      Write output to FILE [stdout]\n");
    fprintf(fp, "  -q, --min-BQ INT\n"
                "               Filter bases with base quality smaller than INT [0]\n");
    fprintf(fp, "  -Q, --min-MQ INT\n"
                "               Filter alignments with mapping quality smaller than INT [0]\n");
    fprintf(fp, "  -J           Include reads with deletions in depth computation\n");
    fprintf(fp, "  -s           Do not count overlapping reads within a template\n");
    sam_global_opt_help(fp, "-.--.@-.");
    exit(exit_status);
}

int main_depth(int argc, char *argv[])
{
    int nfiles, i;
    samFile **fp;
    sam_hdr_t **header;
    int c, has_index_file = 0;
    char *file_list = NULL, **fn = NULL;
    char *out_file = NULL;
    depth_opt opt = {
        .flag = BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL,
        .incl_flag = 0,
        .require_flag = 0,
        .min_qual = 0,
        .min_mqual = 0,
        .skip_del = 1,
        .header = 0,
        .min_len = 0,
        .out = stdout,
        .all_pos = 0,
        .remove_overlaps = 0,
        .reg = NULL,
        .bed = NULL,
    };

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        {"min-MQ",        required_argument, NULL, 'Q'},
        {"min-mq",        required_argument, NULL, 'Q'},
        {"min-BQ",        required_argument, NULL, 'q'},
        {"min-bq",        required_argument, NULL, 'q'},
        {"excl-flags",    required_argument, NULL, 'G'},
        {"incl-flags",    required_argument, NULL, 1},
        {"require-flags", required_argument, NULL, 2},
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '@'),
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "@:q:Q:JHd:m:l:g:G:o:ar:Xf:b:s",
                            lopts, NULL)) >= 0) {
        switch (c) {
        case 'a':
            opt.all_pos++;
            break;

        case 'b':
            opt.bed = bed_read(optarg);
            if (!opt.bed) {
                print_error_errno("depth", "Could not read file \"%s\"",
                                  optarg);
                return 1;
            }
            break;

        case 'f':
            file_list = optarg;
            break;

        case 'd':
        case 'm':
            // depth limit - now ignored
            break;

        case 'g':
            opt.flag &= ~bam_str2flag(optarg);
            break;
        case 'G': // reject if any set
            opt.flag |= bam_str2flag(optarg);
            break;
        case 1: // reject unless at least one set (0 means ignore option)
            opt.incl_flag |= bam_str2flag(optarg);
            break;
        case 2: // reject unless all set
            opt.require_flag |= bam_str2flag(optarg);
            break;

        case 'l':
            opt.min_len = atoi(optarg);
            break;

        case 'H':
            opt.header = 1;
            break;

        case 'q':
            opt.min_qual = atoi(optarg);
            break;
        case 'Q':
            opt.min_mqual = atoi(optarg);
            break;

        case 'J':
            opt.skip_del = 0;
            break;

        case 'o':
            if (opt.out != stdout)
                break;
            opt.out = fopen(out_file = optarg, "w");
            if (!opt.out) {
                print_error_errno("depth", "Cannot open \"%s\" for writing.",
                                  optarg);
                return EXIT_FAILURE;
            }
            break;

        case 'r':
            opt.reg = optarg;
            break;

        case 's':
            opt.remove_overlaps = 1;
            break;

        case 'X':
            has_index_file = 1;
            break;

        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    if (argc < optind+1 && !file_list) {
        if (argc == optind)
            usage_exit(stdout, EXIT_SUCCESS);
        else
            usage_exit(stderr, EXIT_FAILURE);
    }

    if (file_list) {
        if (has_index_file) {
            print_error("depth", "The -f option cannot be combined with -X");
            return 1;
        }
        if (read_file_list(file_list, &nfiles, &fn))
            return 1;
        argv = fn;
        argc = nfiles;
        optind = 0;
    } else {
        nfiles = argc - optind;
    }

    if (has_index_file) {
        if (nfiles%1) {
            print_error("depth", "-X needs one index specified per bam file");
            return 1;
        }
        nfiles /= 2;
    }
    fp = malloc(nfiles * sizeof(*fp));
    header = malloc(nfiles * sizeof(*header));
    if (!fp || !header) {
        print_error_errno("depth", "Out of memory");
        return 1;
    }

    hts_itr_t **itr = NULL;
    if (opt.reg) {
        itr = calloc(nfiles, sizeof(*itr));
        if (!itr)
            return 1;
    }

    for (i = 0; i < nfiles; i++, optind++) {
        fp[i] = sam_open_format(argv[optind], "r", &ga.in);
        if (fp[i] == NULL) {
            print_error_errno("depth",
                              "Cannot open input file \"%s\"", argv[optind]);
            return 1;
        }

        if (ga.nthreads > 0)
            hts_set_threads(fp[i], ga.nthreads);

        if (hts_set_opt(fp[i], CRAM_OPT_REQUIRED_FIELDS,
                        SAM_FLAG | SAM_RNAME | SAM_POS | SAM_CIGAR
                        | (opt.remove_overlaps ? SAM_QNAME|SAM_RNEXT|SAM_PNEXT
                                               : 0)
                        | (opt.min_mqual       ? SAM_MAPQ  : 0)
                        | (opt.min_len         ? SAM_SEQ   : 0)
                        | (opt.min_qual        ? SAM_QUAL  : 0))) {
            fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
            return 1;
        }

        if (hts_set_opt(fp[i], CRAM_OPT_DECODE_MD, 0)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            return 1;
        }

        // FIXME: what if headers differ?
        header[i] = sam_hdr_read(fp[i]);
        if (header == NULL) {
            fprintf(stderr, "Failed to read header for \"%s\"\n",
                    argv[optind]);
            return 1;
        }

        if (opt.reg) {
            hts_idx_t *idx = has_index_file
                ? sam_index_load2(fp[i], argv[optind], argv[optind+nfiles])
                : sam_index_load(fp[i], argv[optind]);
            if (!idx) {
                print_error("depth", "cannot load index for \"%s\"",
                            argv[optind]);
                return 1;
            }
            if (!(itr[i] = sam_itr_querys(idx, header[i], opt.reg))) {
                print_error("depth", "cannot parse region \"%s\"", opt.reg);
                return 1;
            }
            hts_idx_destroy(idx);
        }
    }

    int ret = fastdepth_core(&opt, nfiles, &argv[argc-nfiles], fp, itr, header)
        ? 1 : 0;

    for (i = 0; i < nfiles; i++) {
        sam_hdr_destroy(header[i]);
        sam_close(fp[i]);
        if (itr && itr[i])
            hts_itr_destroy(itr[i]);
    }
    free(header);
    free(fp);
    free(itr);
    if (file_list) {
        for (i=0; i<nfiles; i++)
            free(fn[i]);
        free(fn);
    }
    if (opt.bed)
        bed_destroy(opt.bed);
    sam_global_args_free(&ga);
    if (opt.out != stdout) {
        if (fclose(opt.out) != 0 && ret == 0) {
            print_error_errno("depth", "error on closing \"%s\"", out_file);
            ret = 1;
        }
    }

    return ret;
}

#ifdef _MAIN_BAM2DEPTH
int main(int argc, char *argv[])
{
    return main_depth(argc, argv);
}
#endif
