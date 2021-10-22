/*  bam_aux.c -- remaining aux field handling.

    Copyright (C) 2008-2010, 2013, 2015, 2019 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.

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

#include <ctype.h>
#include "htslib/sam.h"

static inline int bam_aux_type2size(int x)
{
    if (x == 'C' || x == 'c' || x == 'A') return 1;
    else if (x == 'S' || x == 's') return 2;
    else if (x == 'I' || x == 'i' || x == 'f' || x == 'F') return 4;
    else return 0;
}

#define __skip_tag(s) do { \
        int type = toupper(*(s)); \
        ++(s); \
        if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
        else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
        else (s) += bam_aux_type2size(type); \
    } while(0)


int bam_aux_drop_other(bam1_t *b, uint8_t *s)
{
    if (s) {
        uint8_t *p, *aux;
        aux = bam_get_aux(b);
        p = s - 2;
        __skip_tag(s);
        memmove(aux, p, s - p);
        b->l_data -= bam_get_l_aux(b) - (s - p);
    } else {
        b->l_data -= bam_get_l_aux(b);
    }
    return 0;
}
