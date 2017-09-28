/*  bam_lpileup.c -- lplbuf routines.

    Copyright (C) 2008, 2009, 2013 Genome Research Ltd.

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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "bam_plbuf.h"
#include "bam_lpileup.h"
#include "samtools.h"
#include <htslib/ksort.h>

#define TV_GAP 2

typedef struct __freenode_t {
    uint32_t level:28, cnt:4;
    struct __freenode_t *next;
} freenode_t, *freenode_p;

#define freenode_lt(a,b) ((a)->cnt < (b)->cnt || ((a)->cnt == (b)->cnt && (a)->level < (b)->level))
KSORT_INIT(node, freenode_p, freenode_lt)

/* Memory pool, similar to the one in bam_pileup.c */
typedef struct {
    int cnt, n, max;
    freenode_t **buf;
} mempool_t;

static mempool_t *mp_init(void)
{
    return (mempool_t*)calloc(1, sizeof(mempool_t));
}
static void mp_destroy(mempool_t *mp)
{
    int k;
    for (k = 0; k < mp->n; ++k) free(mp->buf[k]);
    free(mp->buf); free(mp);
}
static inline freenode_t *mp_alloc(mempool_t *mp)
{
    ++mp->cnt;
    if (mp->n == 0) return (freenode_t*)calloc(1, sizeof(freenode_t));
    else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, freenode_t *p)
{
    --mp->cnt; p->next = 0; p->cnt = TV_GAP;
    if (mp->n == mp->max) {
        mp->max = mp->max? mp->max<<1 : 256;
        mp->buf = (freenode_t**)realloc(mp->buf, sizeof(freenode_t*) * mp->max);
    }
    mp->buf[mp->n++] = p;
}

/* core part */
struct __bam_lplbuf_t {
    int max, n_cur, n_pre;
    int max_level, *cur_level, *pre_level;
    mempool_t *mp;
    freenode_t **aux, *head, *tail;
    int n_nodes, m_aux;
    bam_pileup_f func;
    void *user_data;
    bam_plbuf_t *plbuf;
};

void bam_lplbuf_reset(bam_lplbuf_t *buf)
{
    freenode_t *p, *q;
    bam_plbuf_reset(buf->plbuf);
    for (p = buf->head; p->next;) {
        q = p->next;
        mp_free(buf->mp, p);
        p = q;
    }
    buf->head = buf->tail;
    buf->max_level = 0;
    buf->n_cur = buf->n_pre = 0;
    buf->n_nodes = 0;
}

static int tview_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    bam_lplbuf_t *tv = (bam_lplbuf_t*)data;
    freenode_t *p;
    int i, l, max_level;
    // allocate memory if necessary
    if (tv->max < n) { // enlarge
        tv->max = n;
        kroundup32(tv->max);
        tv->cur_level = (int*)realloc(tv->cur_level, sizeof(int) * tv->max);
        tv->pre_level = (int*)realloc(tv->pre_level, sizeof(int) * tv->max);
    }
    tv->n_cur = n;
    // update cnt
    for (p = tv->head; p->next; p = p->next)
        if (p->cnt > 0) --p->cnt;
    // calculate cur_level[]
    max_level = 0;
    for (i = l = 0; i < n; ++i) {
        const bam_pileup1_t *p = pl + i;
        if (p->is_head) {
            if (tv->head->next && tv->head->cnt == 0) { // then take a free slot
                freenode_t *p = tv->head->next;
                tv->cur_level[i] = tv->head->level;
                mp_free(tv->mp, tv->head);
                tv->head = p;
                --tv->n_nodes;
            } else tv->cur_level[i] = ++tv->max_level;
        } else {
            tv->cur_level[i] = tv->pre_level[l++];
            if (p->is_tail) { // then return a free slot
                tv->tail->level = tv->cur_level[i];
                tv->tail->next = mp_alloc(tv->mp);
                tv->tail = tv->tail->next;
                ++tv->n_nodes;
            }
        }
        if (tv->cur_level[i] > max_level) max_level = tv->cur_level[i];
        ((bam_pileup1_t*)p)->level = tv->cur_level[i];
    }
    assert(l == tv->n_pre);
    tv->func(tid, pos, n, pl, tv->user_data);
    // sort the linked list
    if (tv->n_nodes) {
        freenode_t *q;
        if (tv->n_nodes + 1 > tv->m_aux) { // enlarge
            tv->m_aux = tv->n_nodes + 1;
            kroundup32(tv->m_aux);
            tv->aux = (freenode_t**)realloc(tv->aux, sizeof(freenode_t*) * tv->m_aux);
        }
        for (p = tv->head, i = l = 0; p->next;) {
            if (p->level > max_level) { // then discard this entry
                q = p->next;
                mp_free(tv->mp, p);
                p = q;
            } else {
                tv->aux[i++] = p;
                p = p->next;
            }
        }
        tv->aux[i] = tv->tail; // add a proper tail for the loop below
        tv->n_nodes = i;
        if (tv->n_nodes) {
            ks_introsort(node, tv->n_nodes, tv->aux);
            for (i = 0; i < tv->n_nodes; ++i) tv->aux[i]->next = tv->aux[i+1];
            tv->head = tv->aux[0];
        } else tv->head = tv->tail;
    }
    // clean up
    tv->max_level = max_level;
    memcpy(tv->pre_level, tv->cur_level, tv->n_cur * 4);
    // squeeze out terminated levels
    for (i = l = 0; i < n; ++i) {
        const bam_pileup1_t *p = pl + i;
        if (!p->is_tail)
            tv->pre_level[l++] = tv->pre_level[i];
    }
    tv->n_pre = l;
/*
    fprintf(stderr, "%d\t", pos+1);
    for (i = 0; i < n; ++i) {
        const bam_pileup1_t *p = pl + i;
        if (p->is_head) fprintf(stderr, "^");
        if (p->is_tail) fprintf(stderr, "$");
        fprintf(stderr, "%d,", p->level);
    }
    fprintf(stderr, "\n");
*/
    return 0;
}

bam_lplbuf_t *bam_lplbuf_init(bam_pileup_f func, void *data)
{
    bam_lplbuf_t *tv;
    tv = (bam_lplbuf_t*)calloc(1, sizeof(bam_lplbuf_t));
    tv->mp = mp_init();
    tv->head = tv->tail = mp_alloc(tv->mp);
    tv->func = func;
    tv->user_data = data;
    tv->plbuf = bam_plbuf_init(tview_func, tv);
    return (bam_lplbuf_t*)tv;
}

void bam_lplbuf_destroy(bam_lplbuf_t *tv)
{
    freenode_t *p, *q;
    free(tv->cur_level); free(tv->pre_level);
    bam_plbuf_destroy(tv->plbuf);
    free(tv->aux);
    for (p = tv->head; p->next;) {
        q = p->next;
        mp_free(tv->mp, p); p = q;
    }
    mp_free(tv->mp, p);
    assert(tv->mp->cnt == 0);
    mp_destroy(tv->mp);
    free(tv);
}

int bam_lplbuf_push(const bam1_t *b, bam_lplbuf_t *tv)
{
    return bam_plbuf_push(b, tv->plbuf);
}
