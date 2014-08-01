/*  bam_plbuf.h -- plbuf routines (declarations copied from bam.h).

    Copyright (C) 2008, 2013 Genome Research Ltd.

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

#ifndef BAM_PLBUF_H
#define BAM_PLBUF_H

#include <htslib/sam.h>

#ifndef BAM_PILEUP_F_DEFINED
#define BAM_PILEUP_F_DEFINED
typedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
#endif //BAM_PILEUP_F_DEFINED

typedef struct {
    bam_plp_t iter;
    bam_pileup_f func;
    void *data;
} bam_plbuf_t;

#ifdef __cplusplus
extern "C" {
#endif
    void bam_plbuf_reset(bam_plbuf_t *buf);

    bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data);

    void bam_plbuf_destroy(bam_plbuf_t *buf);

    int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf);
#ifdef __cplusplus
}
#endif

#endif // BAM_PLBUF_H
