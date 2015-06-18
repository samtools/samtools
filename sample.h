/*  sample.h -- group data by sample.

    Copyright (C) 2010 Broad Institute.

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

#ifndef BAM_SAMPLE_H
#define BAM_SAMPLE_H

#include "htslib/kstring.h"

typedef struct {
    int n, m;
    char **smpl;
    void *rg2smid, *sm2id;
} bam_sample_t;

bam_sample_t *bam_smpl_init(void);
int bam_smpl_add(bam_sample_t *sm, const char *abs, const char *txt);
int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str);
void bam_smpl_destroy(bam_sample_t *sm);

#endif
