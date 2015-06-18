/*  bam_rmdupse.c -- duplicate read detection for unpaired reads.

    Copyright (C) 2009 Genome Research Ltd.
    Portions copyright (C) 2009 Broad Institute.

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

#include <math.h>
#include "sam.h"
#include "htslib/khash.h"
#include "htslib/klist.h"

#define QUEUE_CLEAR_SIZE 0x100000
#define MAX_POS 0x7fffffff

typedef struct {
    int endpos;
    uint32_t score:31, discarded:1;
    bam1_t *b;
} elem_t, *elem_p;
#define __free_elem(p) bam_destroy1((p)->data.b)
KLIST_INIT(q, elem_t, __free_elem)
typedef klist_t(q) queue_t;

KHASH_MAP_INIT_INT(best, elem_p)
typedef khash_t(best) besthash_t;

typedef struct {
    uint64_t n_checked, n_removed;
    besthash_t *left, *rght;
} lib_aux_t;
KHASH_MAP_INIT_STR(lib, lib_aux_t)

static lib_aux_t *get_aux(khash_t(lib) *aux, const char *lib)
{
    khint_t k = kh_get(lib, aux, lib);
    if (k == kh_end(aux)) {
        int ret;
        char *p = strdup(lib);
        lib_aux_t *q;
        k = kh_put(lib, aux, p, &ret);
        q = &kh_val(aux, k);
        q->left = kh_init(best);
        q->rght = kh_init(best);
        q->n_checked = q->n_removed = 0;
        return q;
    } else return &kh_val(aux, k);
}

static inline int sum_qual(const bam1_t *b)
{
    int i, q;
    uint8_t *qual = bam1_qual(b);
    for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
    return q;
}

static inline elem_t *push_queue(queue_t *queue, const bam1_t *b, int endpos, int score)
{
    elem_t *p = kl_pushp(q, queue);
    p->discarded = 0;
    p->endpos = endpos; p->score = score;
    if (p->b == 0) p->b = bam_init1();
    bam_copy1(p->b, b);
    return p;
}

static void clear_besthash(besthash_t *h, int32_t pos)
{
    khint_t k;
    for (k = kh_begin(h); k != kh_end(h); ++k)
        if (kh_exist(h, k) && kh_val(h, k)->endpos <= pos)
            kh_del(best, h, k);
}

static void dump_alignment(samfile_t *out, queue_t *queue, int32_t pos, khash_t(lib) *h)
{
    if (queue->size > QUEUE_CLEAR_SIZE || pos == MAX_POS) {
        khint_t k;
        while (1) {
            elem_t *q;
            if (queue->head == queue->tail) break;
            q = &kl_val(queue->head);
            if (q->discarded) {
                q->b->data_len = 0;
                kl_shift(q, queue, 0);
                continue;
            }
            if ((q->b->core.flag&BAM_FREVERSE) && q->endpos > pos) break;
            samwrite(out, q->b);
            q->b->data_len = 0;
            kl_shift(q, queue, 0);
        }
        for (k = kh_begin(h); k != kh_end(h); ++k) {
            if (kh_exist(h, k)) {
                clear_besthash(kh_val(h, k).left, pos);
                clear_besthash(kh_val(h, k).rght, pos);
            }
        }
    }
}

void bam_rmdupse_core(samfile_t *in, samfile_t *out, int force_se)
{
    bam1_t *b;
    queue_t *queue;
    khint_t k;
    int last_tid = -2;
    khash_t(lib) *aux;

    aux = kh_init(lib);
    b = bam_init1();
    queue = kl_init(q);
    while (samread(in, b) >= 0) {
        bam1_core_t *c = &b->core;
        int endpos = bam_calend(c, bam1_cigar(b));
        int score = sum_qual(b);

        if (last_tid != c->tid) {
            if (last_tid >= 0) dump_alignment(out, queue, MAX_POS, aux);
            last_tid = c->tid;
        } else dump_alignment(out, queue, c->pos, aux);
        if ((c->flag&BAM_FUNMAP) || ((c->flag&BAM_FPAIRED) && !force_se)) {
            push_queue(queue, b, endpos, score);
        } else {
            const char *lib;
            lib_aux_t *q;
            besthash_t *h;
            uint32_t key;
            int ret;
            lib = bam_get_library(in->header, b);
            q = lib? get_aux(aux, lib) : get_aux(aux, "\t");
            ++q->n_checked;
            h = (c->flag&BAM_FREVERSE)? q->rght : q->left;
            key = (c->flag&BAM_FREVERSE)? endpos : c->pos;
            k = kh_put(best, h, key, &ret);
            if (ret == 0) { // in the hash table
                elem_t *p = kh_val(h, k);
                ++q->n_removed;
                if (p->score < score) {
                    if (c->flag&BAM_FREVERSE) { // mark "discarded" and push the queue
                        p->discarded = 1;
                        kh_val(h, k) = push_queue(queue, b, endpos, score);
                    } else { // replace
                        p->score = score; p->endpos = endpos;
                        bam_copy1(p->b, b);
                    }
                } // otherwise, discard the alignment
            } else kh_val(h, k) = push_queue(queue, b, endpos, score);
        }
    }
    dump_alignment(out, queue, MAX_POS, aux);

    for (k = kh_begin(aux); k != kh_end(aux); ++k) {
        if (kh_exist(aux, k)) {
            lib_aux_t *q = &kh_val(aux, k);
            fprintf(stderr, "[bam_rmdupse_core] %lld / %lld = %.4lf in library '%s'\n", (long long)q->n_removed,
                    (long long)q->n_checked, (double)q->n_removed/q->n_checked, kh_key(aux, k));
            kh_destroy(best, q->left); kh_destroy(best, q->rght);
            free((char*)kh_key(aux, k));
        }
    }
    kh_destroy(lib, aux);
    bam_destroy1(b);
    kl_destroy(q, queue);
}
