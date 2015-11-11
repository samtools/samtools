/*  bam2bcf.h -- variant calling.

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

#ifndef BAM2BCF_H
#define BAM2BCF_H

#include <stdint.h>
#include <htslib/vcf.h>
#include "errmod.h"

/**
 *  A simplified version of Mann-Whitney U-test is calculated
 *  by default (no CDF) because it is faster and seems to work
 *  better in machine learning filtering. When enabled by setting
 *  CDF_MWU_TESTS, additional annotations will appear on mpileup's
 *  output (RPB2 in addition to RPB, etc.).
 */
#ifndef CDF_MWU_TESTS
#define CDF_MWU_TESTS 0
#endif

#define B2B_INDEL_NULL 10000

#define B2B_FMT_DP      (1<<0)
#define B2B_FMT_SP      (1<<1)
#define B2B_FMT_DV      (1<<2)
#define B2B_FMT_DP4     (1<<3)
#define B2B_FMT_DPR     (1<<4)
#define B2B_INFO_DPR    (1<<5)
#define B2B_FMT_AD      (1<<6)
#define B2B_FMT_ADF     (1<<7)
#define B2B_FMT_ADR     (1<<8)
#define B2B_INFO_AD     (1<<9)
#define B2B_INFO_ADF    (1<<10)
#define B2B_INFO_ADR    (1<<11)

#define B2B_MAX_ALLELES 5

typedef struct __bcf_callaux_t {
    int capQ, min_baseQ;
    int openQ, extQ, tandemQ; // for indels
    uint32_t min_support, max_support; // for collecting indel candidates
    double min_frac; // for collecting indel candidates
    float max_frac; // for collecting indel candidates
    int per_sample_flt; // indel filtering strategy
    int *ref_pos, *alt_pos, npos, *ref_mq, *alt_mq, *ref_bq, *alt_bq, *fwd_mqs, *rev_mqs, nqual; // for bias tests
    // for internal uses
    int max_bases;
    int indel_types[4];     // indel lengths
    int maxins, indelreg;
    int read_len;
    char *inscns;
    uint16_t *bases;        // 5bit: unused, 6:quality, 1:is_rev, 4:2-bit base or indel allele (index to bcf_callaux_t.indel_types)
    errmod_t *e;
    void *rghash;
} bcf_callaux_t;

typedef struct {
    uint32_t ori_depth;
    unsigned int mq0;
    int32_t *ADF, *ADR;
    float qsum[4];
    // The fields are:
    //      depth fwd   .. ref (0) and non-ref (2)
    //      depth rev   .. ref (1) and non-ref (3)
    //      baseQ       .. ref (4) and non-ref (6)
    //      baseQ^2     .. ref (5) and non-ref (7)
    //      mapQ        .. ref (8) and non-ref (10)
    //      mapQ^2      .. ref (9) and non-ref (11)
    //      minDist     .. ref (12) and non-ref (14)
    //      minDist^2   .. ref (13) and non-ref (15)
    // Note that this probably needs a more thorough fix: int types in
    // bcf_call_t do overflow with high-coverage data, such as exomes, and
    // BCFv2 supports only floats which may not suffice.
    double anno[16];
    float p[25];        // phred-scaled likelihood of each genotype
} bcf_callret1_t;

typedef struct {
    int tid, pos;
    bcf_hdr_t *bcf_hdr;
    int a[5]; // alleles: ref, alt, alt2, alt3
    float qsum[5];  // for the QS tag
    int n, n_alleles, shift, ori_ref, unseen;
    int n_supp; // number of supporting non-reference reads
    double anno[16];
    unsigned int depth, ori_depth, mq0;
    int32_t *PL, *DP4, *ADR, *ADF;
    uint8_t *fmt_arr;
    float vdb; // variant distance bias
    float mwu_pos, mwu_mq, mwu_bq, mwu_mqs;
#if CDF_MWU_TESTS
    float mwu_pos_cdf, mwu_mq_cdf, mwu_bq_cdf, mwu_mqs_cdf;
#endif
    float seg_bias;
    kstring_t tmp;
} bcf_call_t;

#ifdef __cplusplus
extern "C" {
#endif

    bcf_callaux_t *bcf_call_init(double theta, int min_baseQ);
    void bcf_call_destroy(bcf_callaux_t *bca);
    int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r);
    int bcf_call_combine(int n, const bcf_callret1_t *calls, bcf_callaux_t *bca, int ref_base /*4-bit*/, bcf_call_t *call);
    int bcf_call2bcf(bcf_call_t *bc, bcf1_t *b, bcf_callret1_t *bcr, int fmt_flag,
                     const bcf_callaux_t *bca, const char *ref);
    int bcf_call_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref,
                          const void *rghash);
    void bcf_callaux_clean(bcf_callaux_t *bca, bcf_call_t *call);

#ifdef __cplusplus
}
#endif

#endif
