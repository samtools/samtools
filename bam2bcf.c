/*  bam2bcf.c -- variant calling.

    Copyright (C) 2010-2012 Broad Institute.
    Copyright (C) 2012-2014 Genome Research Ltd.

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

#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <float.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/kfunc.h>
#include "bam2bcf.h"
#include "errmod.h"

extern  void ks_introsort_uint32_t(size_t n, uint32_t a[]);
extern const char bam_nt16_nt4_table[];

#define CALL_DEFTHETA 0.83
#define DEF_MAPQ 20

#define CAP_DIST 25

bcf_callaux_t *bcf_call_init(double theta, int min_baseQ)
{
    bcf_callaux_t *bca;
    if (theta <= 0.) theta = CALL_DEFTHETA;
    bca = calloc(1, sizeof(bcf_callaux_t));
    bca->capQ = 60;
    bca->openQ = 40; bca->extQ = 20; bca->tandemQ = 100;
    bca->min_baseQ = min_baseQ;
    bca->e = errmod_init(1. - theta);
    bca->min_frac = 0.002;
    bca->min_support = 1;
    bca->per_sample_flt = 0;
    bca->npos = 100;
    bca->ref_pos = malloc(bca->npos*sizeof(int));
    bca->alt_pos = malloc(bca->npos*sizeof(int));
    bca->nqual = 60;
    bca->ref_mq  = malloc(bca->nqual*sizeof(int));
    bca->alt_mq  = malloc(bca->nqual*sizeof(int));
    bca->ref_bq  = malloc(bca->nqual*sizeof(int));
    bca->alt_bq  = malloc(bca->nqual*sizeof(int));
    bca->fwd_mqs = malloc(bca->nqual*sizeof(int));
    bca->rev_mqs = malloc(bca->nqual*sizeof(int));
    return bca;
}

void bcf_call_destroy(bcf_callaux_t *bca)
{
    if (bca == 0) return;
    errmod_destroy(bca->e);
    if (bca->npos) { free(bca->ref_pos); free(bca->alt_pos); bca->npos = 0; }
    free(bca->ref_mq); free(bca->alt_mq); free(bca->ref_bq); free(bca->alt_bq);
    free(bca->fwd_mqs); free(bca->rev_mqs);
    bca->nqual = 0;
    free(bca->bases); free(bca->inscns); free(bca);
}

// position in the sequence with respect to the aligned part of the read
static int get_position(const bam_pileup1_t *p, int *len)
{
    int icig, n_tot_bases = 0, iread = 0, edist = p->qpos + 1;
    for (icig=0; icig<p->b->core.n_cigar; icig++)
    {
        int cig  = bam_get_cigar(p->b)[icig] & BAM_CIGAR_MASK;
        int ncig = bam_get_cigar(p->b)[icig] >> BAM_CIGAR_SHIFT;
        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            n_tot_bases += ncig;
            iread += ncig;
            continue;
        }
        if ( cig==BAM_CINS )
        {
            n_tot_bases += ncig;
            iread += ncig;
            continue;
        }
        if ( cig==BAM_CSOFT_CLIP )
        {
            iread += ncig;
            if ( iread<=p->qpos ) edist -= ncig;
            continue;
        }
        if ( cig==BAM_CDEL ) continue;
        if ( cig==BAM_CHARD_CLIP ) continue;
        if ( cig==BAM_CPAD ) continue;
        if ( cig==BAM_CREF_SKIP ) continue;
        fprintf(stderr,"todo: cigar %d\n", cig);
        assert(0);
    }
    *len = n_tot_bases;
    return edist;
}

/*
    Notes:
    - Called from bam_plcmd.c by mpileup. Amongst other things, sets the bcf_callret1_t.qsum frequencies
        which are carried over via bcf_call_combine and bcf_call2bcf to the output BCF as the QS annotation.
        Later it's used for multiallelic calling by bcftools -m
    - ref_base is the 4-bit representation of the reference base. It is negative if we are looking at an indel.
 */
/*
 * This function is called once for each sample.
 * _n is number of pilesups pl contributing reads to this sample
 * pl is pointer to array of _n pileups (one pileup per read)
 * ref_base is the 4-bit representation of the reference base. It is negative if we are looking at an indel.
 * bca is the settings to perform calls across all samples
 * r is the returned value of the call
 */
int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r)
{
    int i, n, ref4, is_indel, ori_depth = 0;

    // clean from previous run
    r->ori_depth = 0;
    r->mq0 = 0;
    memset(r->qsum,0,sizeof(float)*4);
    memset(r->anno,0,sizeof(double)*16);
    memset(r->p,0,sizeof(float)*25);

    if (ref_base >= 0) {
        ref4 = bam_nt16_nt4_table[ref_base];
        is_indel = 0;
    } else ref4 = 4, is_indel = 1;
    if (_n == 0) return -1;
    // enlarge the bases array if necessary
    if (bca->max_bases < _n) {
        bca->max_bases = _n;
        kroundup32(bca->max_bases);
        bca->bases = (uint16_t*)realloc(bca->bases, 2 * bca->max_bases);
    }
    // fill the bases array
    for (i = n = 0; i < _n; ++i) {
        const bam_pileup1_t *p = pl + i;
        int q, b, mapQ, baseQ, is_diff, min_dist, seqQ;
        // set base
        if (p->is_del || p->is_refskip || (p->b->core.flag&BAM_FUNMAP)) continue;
        ++ori_depth;
        mapQ  = p->b->core.qual < 255? p->b->core.qual : DEF_MAPQ; // special case for mapQ==255
        if ( !mapQ ) r->mq0++;
        baseQ = q = is_indel? p->aux&0xff : (int)bam_get_qual(p->b)[p->qpos]; // base/indel quality
        seqQ = is_indel? (p->aux>>8&0xff) : 99;
        if (q < bca->min_baseQ) continue;
        if (q > seqQ) q = seqQ;
        mapQ = mapQ < bca->capQ? mapQ : bca->capQ;
        if (q > mapQ) q = mapQ;
        if (q > 63) q = 63;
        if (q < 4) q = 4;       // MQ=0 reads count as BQ=4
        if (!is_indel) {
            b = bam_seqi(bam_get_seq(p->b), p->qpos); // base
            b = bam_nt16_nt4_table[b? b : ref_base]; // b is the 2-bit base
            is_diff = (ref4 < 4 && b == ref4)? 0 : 1;
        } else {
            b = p->aux>>16&0x3f;
            is_diff = (b != 0);
        }
        bca->bases[n++] = q<<5 | (int)bam_is_rev(p->b)<<4 | b;
        // collect annotations
        if (b < 4)
        {
            r->qsum[b] += q;
            if ( r->DPR ) r->DPR[b]++;
        }
        ++r->anno[0<<2|is_diff<<1|bam_is_rev(p->b)];
        min_dist = p->b->core.l_qseq - 1 - p->qpos;
        if (min_dist > p->qpos) min_dist = p->qpos;
        if (min_dist > CAP_DIST) min_dist = CAP_DIST;
        r->anno[1<<2|is_diff<<1|0] += baseQ;
        r->anno[1<<2|is_diff<<1|1] += baseQ * baseQ;
        r->anno[2<<2|is_diff<<1|0] += mapQ;
        r->anno[2<<2|is_diff<<1|1] += mapQ * mapQ;
        r->anno[3<<2|is_diff<<1|0] += min_dist;
        r->anno[3<<2|is_diff<<1|1] += min_dist * min_dist;

        // collect for bias tests
        if ( baseQ > 59 ) baseQ = 59;
        if ( mapQ > 59 ) mapQ = 59;
        int len, pos = get_position(p, &len);
        int epos = (double)pos/(len+1) * bca->npos;
        int ibq  = baseQ/60. * bca->nqual;
        int imq  = mapQ/60. * bca->nqual;
        if ( bam_is_rev(p->b) ) bca->rev_mqs[imq]++;
        else bca->fwd_mqs[imq]++;
        if ( bam_seqi(bam_get_seq(p->b),p->qpos) == ref_base )
        {
            bca->ref_pos[epos]++;
            bca->ref_bq[ibq]++;
            bca->ref_mq[imq]++;
        }
        else
        {
            bca->alt_pos[epos]++;
            bca->alt_bq[ibq]++;
            bca->alt_mq[imq]++;
        }
    }
    r->ori_depth = ori_depth;
    // glfgen
    errmod_cal(bca->e, n, 5, bca->bases, r->p); // calculate PL of each genotype
    return n;
}
