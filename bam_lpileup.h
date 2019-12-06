/*  bam_lpileup.h -- lplbuf routines (declarations copied from bam.h).

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

#ifndef BAM_LPILEUP_H
#define BAM_LPILEUP_H


#include <htslib/sam.h>

struct __bam_lplbuf_t;
typedef struct __bam_lplbuf_t bam_lplbuf_t;

#ifndef BAM_PILEUP_F_DEFINED
#define BAM_PILEUP_F_DEFINED
typedef int (*bam_pileup_f)(uint32_t tid, hts_pos_t pos, int n, const bam_pileup1_t *pl, void *data);
#endif //BAM_PILEUP_F_DEFINED


#ifdef __cplusplus
extern "C" {
#endif
    void bam_lplbuf_reset(bam_lplbuf_t *buf);

    /*! @abstract  bam_plbuf_init() equivalent with level calculated. */
    bam_lplbuf_t *bam_lplbuf_init(bam_pileup_f func, void *data);

    /*! @abstract  bam_plbuf_destroy() equivalent with level calculated. */
    void bam_lplbuf_destroy(bam_lplbuf_t *tv);

    /*! @abstract  bam_plbuf_push() equivalent with level calculated. */
    int bam_lplbuf_push(const bam1_t *b, bam_lplbuf_t *buf);
#ifdef __cplusplus
}
#endif

#endif // BAM_LPILEUP_H
