/*  bam2bcf.c -- variant calling.  Used for tview consensus

    Copyright (C) 2010-2012 Broad Institute.
    Copyright (C) 2012-2015, 2021 Genome Research Ltd.

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

#include <stdint.h>
#include <assert.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include "bam2bcf.h"

#define CALL_DEFTHETA 0.83
#define DEF_MAPQ 20

#define CAP_DIST 25

bcf_callaux_t *bcf_call_init(double theta, int min_baseQ)
{
    bcf_callaux_t *bca;
    if (theta <= 0.) theta = CALL_DEFTHETA;
    bca = calloc(1, sizeof(bcf_callaux_t));
    bca->capQ = 60;
    bca->min_baseQ = min_baseQ;
    bca->e = errmod_init(1. - theta);
    return bca;
}

void bcf_call_destroy(bcf_callaux_t *bca)
{
    if (bca == 0) return;
    errmod_destroy(bca->e);
    free(bca->bases); free(bca);
}

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
    int i, n, is_indel;

    // clean from previous run
    memset(r->qsum,0,sizeof(float)*4);
    memset(r->p,0,sizeof(float)*25);

    if (ref_base >= 0) {
        is_indel = 0;
    } else is_indel = 1;
    if (_n <= 0) return -1;
    // enlarge the bases array if necessary
    if (bca->max_bases < _n) {
        bca->max_bases = _n;
        kroundup32(bca->max_bases);
        bca->bases = (uint16_t*)realloc(bca->bases, 2 * (size_t) bca->max_bases);
    }
    // fill the bases array
    for (i = n = 0; i < _n; ++i) {
        const bam_pileup1_t *p = pl + i;
        int q, b, mapQ, min_dist, seqQ;
        // set base
        if (p->is_del || p->is_refskip || (p->b->core.flag&BAM_FUNMAP)) continue;
        mapQ  = p->b->core.qual < 255? p->b->core.qual : DEF_MAPQ; // special case for mapQ==255
        q = is_indel? p->aux&0xff : (int)bam_get_qual(p->b)[p->qpos]; // base/indel quality
        seqQ = is_indel? (p->aux>>8&0xff) : 99;
        if (q < bca->min_baseQ) continue;
        if (q > seqQ) q = seqQ;
        mapQ = mapQ < bca->capQ? mapQ : bca->capQ;
        if (q > mapQ) q = mapQ;
        if (q > 63) q = 63;
        if (q < 4) q = 4;       // MQ=0 reads count as BQ=4
        if (!is_indel) {
            b = bam_seqi(bam_get_seq(p->b), p->qpos); // base
            b = seq_nt16_int[b? b : ref_base]; // b is the 2-bit base
        } else {
            b = p->aux>>16&0x3f;
        }
        bca->bases[n++] = q<<5 | (int)bam_is_rev(p->b)<<4 | b;
        // collect annotations
        if (b < 4)
            r->qsum[b] += q;

        min_dist = p->b->core.l_qseq - 1 - p->qpos;
        if (min_dist > p->qpos) min_dist = p->qpos;
        if (min_dist > CAP_DIST) min_dist = CAP_DIST;
    }
    // glfgen
    errmod_cal(bca->e, n, 5, bca->bases, r->p); // calculate PL of each genotype
    return n;
}
