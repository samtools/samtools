/*  sample.c -- group data by sample.

    Copyright (C) 2010, 2011 Broad Institute.
    Copyright (C) 2013 Genome Research Ltd.

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

#include <stdlib.h>
#include <string.h>
#include "sample.h"
#include "htslib/khash.h"
KHASH_MAP_INIT_STR(sm, int)

bam_sample_t *bam_smpl_init(void)
{
    bam_sample_t *s;
    s = calloc(1, sizeof(bam_sample_t));
    s->rg2smid = kh_init(sm);
    s->sm2id = kh_init(sm);
    return s;
}

void bam_smpl_destroy(bam_sample_t *sm)
{
    int i;
    khint_t k;
    khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;
    if (sm == 0) return;
    for (i = 0; i < sm->n; ++i) free(sm->smpl[i]);
    free(sm->smpl);
    for (k = kh_begin(rg2smid); k != kh_end(rg2smid); ++k)
        if (kh_exist(rg2smid, k)) free((char*)kh_key(rg2smid, k));
    kh_destroy(sm, sm->rg2smid);
    kh_destroy(sm, sm->sm2id);
    free(sm);
}

static void add_pair(bam_sample_t *sm, khash_t(sm) *sm2id, const char *key, const char *val)
{
    khint_t k_rg, k_sm;
    int ret;
    khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;
    k_rg = kh_get(sm, rg2smid, key);
    if (k_rg != kh_end(rg2smid)) return; // duplicated @RG-ID
    k_rg = kh_put(sm, rg2smid, strdup(key), &ret);
    k_sm = kh_get(sm, sm2id, val);
    if (k_sm == kh_end(sm2id)) { // absent
        if (sm->n == sm->m) {
            sm->m = sm->m? sm->m<<1 : 1;
            sm->smpl = realloc(sm->smpl, sizeof(char*) * sm->m);
        }
        sm->smpl[sm->n] = strdup(val);
        k_sm = kh_put(sm, sm2id, sm->smpl[sm->n], &ret);
        kh_val(sm2id, k_sm) = sm->n++;
    }
    kh_val(rg2smid, k_rg) = kh_val(sm2id, k_sm);
}

int bam_smpl_add(bam_sample_t *sm, const char *fn, const char *txt)
{
    const char *p = txt, *q, *r;
    kstring_t buf, first_sm;
    int n = 0;
    khash_t(sm) *sm2id = (khash_t(sm)*)sm->sm2id;
    if (txt == 0) {
        add_pair(sm, sm2id, fn, fn);
        return 0;
    }
    memset(&buf, 0, sizeof(kstring_t));
    memset(&first_sm, 0, sizeof(kstring_t));
    while ((q = strstr(p, "@RG")) != 0) {
        p = q + 3;
        r = q = 0;
        if ((q = strstr(p, "\tID:")) != 0) q += 4;
        if ((r = strstr(p, "\tSM:")) != 0) r += 4;
        if (r && q) {
            char *u, *v;
            int oq, or;
            for (u = (char*)q; *u && *u != '\t' && *u != '\n'; ++u);
            for (v = (char*)r; *v && *v != '\t' && *v != '\n'; ++v);
            oq = *u; or = *v; *u = *v = '\0';
            buf.l = 0; kputs(fn, &buf); kputc('/', &buf); kputs(q, &buf);
            add_pair(sm, sm2id, buf.s, r);
            if ( !first_sm.s )
                kputs(r,&first_sm);
            *u = oq; *v = or;
        } else break;
        p = q > r? q : r;
        ++n;
    }
    if (n == 0) add_pair(sm, sm2id, fn, fn);
    // If there is only one RG tag present in the header and reads are not annotated, don't refuse to work but
    //  use the tag instead.
    else if ( n==1 && first_sm.s )
        add_pair(sm,sm2id,fn,first_sm.s);
    if ( first_sm.s )
        free(first_sm.s);

//  add_pair(sm, sm2id, fn, fn);
    free(buf.s);
    return 0;
}

int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str)
{
    khint_t k;
    khash_t(sm) *rg2smid = (khash_t(sm)*)sm->rg2smid;
    if (rg) {
        str->l = 0;
        kputs(fn, str); kputc('/', str); kputs(rg, str);
        k = kh_get(sm, rg2smid, str->s);
    } else k = kh_get(sm, rg2smid, fn);
    return k == kh_end(rg2smid)? -1 : kh_val(rg2smid, k);
}
