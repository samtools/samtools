/*  bam2bcf.h -- variant calling.

    Copyright (C) 2010-2012 Broad Institute.
    Copyright (C) 2012-2014, 2019, 2021 Genome Research Ltd.

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
#include <htslib/hts.h>
#include <htslib/vcf.h>

typedef struct __bcf_callaux_t {
    int capQ, min_baseQ;
    // for internal uses
    int max_bases;
    uint16_t *bases;        // 5bit: unused, 6:quality, 1:is_rev, 4:2-bit base or indel allele (index to bcf_callaux_t.indel_types)
    errmod_t *e;
} bcf_callaux_t;

// NB: Only bits used by tview are retained.
typedef struct {
    float qsum[4];
    float p[25];        // phred-scaled likelihood of each genotype
} bcf_callret1_t;

#ifdef __cplusplus
extern "C" {
#endif

    bcf_callaux_t *bcf_call_init(double theta, int min_baseQ);
    void bcf_call_destroy(bcf_callaux_t *bca);
    int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r);

#ifdef __cplusplus
}
#endif

#endif
