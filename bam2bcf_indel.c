/*  bam2bcf_indel.c -- indel caller.

    Copyright (C) 2010, 2011 Broad Institute.
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

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "htslib/sam.h"
#include "bam2bcf.h"
#include "kaln.h"
#include "kprobaln.h"
#include "htslib/khash.h"
KHASH_SET_INIT_STR(rg)

void *bcf_call_add_rg(void *_hash, const char *hdtext, const char *list)
{
    const char *s, *p, *q, *r, *t;
    khash_t(rg) *hash;
    if (list == 0 || hdtext == 0) return _hash;
    if (_hash == 0) _hash = kh_init(rg);
    hash = (khash_t(rg)*)_hash;
    if ((s = strstr(hdtext, "@RG\t")) == 0) return hash;
    do {
        t = strstr(s + 4, "@RG\t"); // the next @RG
        if ((p = strstr(s, "\tID:")) != 0) p += 4;
        if ((q = strstr(s, "\tPL:")) != 0) q += 4;
        if (p && q && (t == 0 || (p < t && q < t))) { // ID and PL are both present
            int lp, lq;
            char *x;
            for (r = p; *r && *r != '\t' && *r != '\n'; ++r) { }
            lp = r - p;
            for (r = q; *r && *r != '\t' && *r != '\n'; ++r) { }
            lq = r - q;
            x = calloc((lp > lq? lp : lq) + 1, 1);
            for (r = q; *r && *r != '\t' && *r != '\n'; ++r) x[r-q] = *r;
            if (strstr(list, x)) { // insert ID to the hash table
                khint_t k;
                int ret;
                for (r = p; *r && *r != '\t' && *r != '\n'; ++r) x[r-p] = *r;
                x[r-p] = 0;
                k = kh_get(rg, hash, x);
                if (k == kh_end(hash)) k = kh_put(rg, hash, x, &ret);
                else free(x);
            } else free(x);
        }
        s = t;
    } while (s);
    return hash;
}

void bcf_call_del_rghash(void *_hash)
{
    khint_t k;
    khash_t(rg) *hash = (khash_t(rg)*)_hash;
    if (hash == 0) return;
    for (k = kh_begin(hash); k < kh_end(hash); ++k)
        if (kh_exist(hash, k))
            free((char*)kh_key(hash, k));
    kh_destroy(rg, hash);
}
