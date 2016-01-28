/*  bam_plbuf.c -- plbuf routines (previously in bam_pileup.c).

    Copyright (C) 2008-2010, 2013 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include "bam_plbuf.h"

/*****************
 * callback APIs *
 *****************/

void bam_plbuf_reset(bam_plbuf_t *buf)
{
    bam_plp_reset(buf->iter);
}

bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)
{
    bam_plbuf_t *buf;
    buf = calloc(1, sizeof(bam_plbuf_t));
    buf->iter = bam_plp_init(0, 0);
    buf->func = func;
    buf->data = data;
    return buf;
}

void bam_plbuf_destroy(bam_plbuf_t *buf)
{
    bam_plp_destroy(buf->iter);
    free(buf);
}

int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf)
{
    int ret, n_plp, tid, pos;
    const bam_pileup1_t *plp;
    ret = bam_plp_push(buf->iter, b);
    if (ret < 0) return ret;
    while ((plp = bam_plp_next(buf->iter, &tid, &pos, &n_plp)) != 0)
        buf->func(tid, pos, n_plp, plp, buf->data);
    return 0;
}
