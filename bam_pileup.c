#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "sam.h"

typedef struct __linkbuf_t {
	bam1_t b;
	uint32_t beg, end;
	struct __linkbuf_t *next;
} lbnode_t;

/* --- BEGIN: Memory pool */

typedef struct {
	int cnt, n, max;
	lbnode_t **buf;
} mempool_t;

static mempool_t *mp_init()
{
	mempool_t *mp;
	mp = (mempool_t*)calloc(1, sizeof(mempool_t));
	return mp;
}
static void mp_destroy(mempool_t *mp)
{
	int k;
	for (k = 0; k < mp->n; ++k) {
		free(mp->buf[k]->b.data);
		free(mp->buf[k]);
	}
	free(mp->buf);
	free(mp);
}
static inline lbnode_t *mp_alloc(mempool_t *mp)
{
	++mp->cnt;
	if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
	else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
	--mp->cnt; p->next = 0; // clear lbnode_t::next here
	if (mp->n == mp->max) {
		mp->max = mp->max? mp->max<<1 : 256;
		mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
	}
	mp->buf[mp->n++] = p;
}

/* --- END: Memory pool */

/* --- BEGIN: Auxiliary functions */

static inline int resolve_cigar(bam_pileup1_t *p, uint32_t pos)
{
	unsigned k;
	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t x = c->pos, y = 0;
	int ret = 1, is_restart = 1;

	if (c->flag&BAM_FUNMAP) return 0; // unmapped read
	assert(x <= pos); // otherwise a bug
	p->qpos = -1; p->indel = 0; p->is_del = p->is_head = p->is_tail = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = bam1_cigar(b)[k] & BAM_CIGAR_MASK; // operation
		int l = bam1_cigar(b)[k] >> BAM_CIGAR_SHIFT; // length
		if (op == BAM_CMATCH) { // NOTE: this assumes the first and the last operation MUST BE a match or a clip
			if (x + l > pos) { // overlap with pos
				p->indel = p->is_del = 0;
				p->qpos = y + (pos - x);
				if (x == pos && is_restart) p->is_head = 1;
				if (x + l - 1 == pos) { // come to the end of a match
					if (k < c->n_cigar - 1) { // there are additional operation(s)
						uint32_t cigar = bam1_cigar(b)[k+1]; // next CIGAR
						int op_next = cigar&BAM_CIGAR_MASK; // next CIGAR operation
						if (op_next == BAM_CDEL) p->indel = -(int32_t)(cigar>>BAM_CIGAR_SHIFT); // del
						else if (op_next == BAM_CINS) p->indel = cigar>>BAM_CIGAR_SHIFT; // ins
						if (op_next == BAM_CDEL || op_next == BAM_CINS) {
							if (k + 2 < c->n_cigar) op_next = bam1_cigar(b)[k+2]&BAM_CIGAR_MASK;
							else p->is_tail = 1;
						}
						if (op_next == BAM_CSOFT_CLIP || op_next == BAM_CREF_SKIP || op_next == BAM_CHARD_CLIP)
							p->is_tail = 1; // tail
					} else p->is_tail = 1; // this is the last operation; set tail
				}
			}
			x += l; y += l;
		} else if (op == BAM_CDEL) { // then set ->is_del
			if (x + l > pos) {
				p->indel = 0; p->is_del = 1;
				p->qpos = y + (pos - x);
			}
			x += l;
		} else if (op == BAM_CREF_SKIP) x += l;
		else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		is_restart = (op == BAM_CREF_SKIP || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP);
		if (x > pos) {
			if (op == BAM_CREF_SKIP) ret = 0; // then do not put it into pileup at all
			break;
		}
	}
	assert(x > pos); // otherwise a bug
	return ret;
}

/* --- END: Auxiliary functions */

/*******************
 * pileup iterator *
 *******************/

struct __bam_plp_t {
	mempool_t *mp;
	lbnode_t *head, *tail, *dummy;
	int32_t tid, pos, max_tid, max_pos;
	int is_eof, flag_mask, max_plp;
	bam_pileup1_t *plp;
};

bam_plp_t bam_plp_init(void)
{
	bam_plp_t iter;
	iter = calloc(1, sizeof(struct __bam_plp_t));
	iter->mp = mp_init();
	iter->head = iter->tail = mp_alloc(iter->mp);
	iter->dummy = mp_alloc(iter->mp);
	iter->max_tid = iter->max_pos = -1;
	iter->flag_mask = BAM_DEF_MASK;
	return iter;
}

const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_n_plp, int *_tid, int *_pos)
{
	*_n_plp = 0;
	if (iter->is_eof && iter->head->next == 0) return 0;
	while (iter->is_eof || iter->max_tid > iter->tid || (iter->max_tid == iter->tid && iter->max_pos > iter->pos)) {
		int n_plp = 0;
		lbnode_t *p, *q;
		// write iter->plp at iter->pos
		iter->dummy->next = iter->head;
		for (p = iter->head, q = iter->dummy; p->next; q = p, p = p->next) {
			if (p->b.core.tid < iter->tid || (p->b.core.tid == iter->tid && p->end <= iter->pos)) { // then remove
				q->next = p->next; mp_free(iter->mp, p); p = q;
			} else if (p->b.core.tid == iter->tid && p->beg <= iter->pos) { // here: p->end > pos; then add to pileup
				if (n_plp == iter->max_plp) { // then double the capacity
					iter->max_plp = iter->max_plp? iter->max_plp<<1 : 256;
					iter->plp = (bam_pileup1_t*)realloc(iter->plp, sizeof(bam_pileup1_t) * iter->max_plp);
				}
				iter->plp[n_plp].b = &p->b;
				if (resolve_cigar(iter->plp + n_plp, iter->pos)) ++n_plp; // skip the read if we are looking at ref-skip
			}
		}
		iter->head = iter->dummy->next; // dummy->next may be changed
		*_n_plp = n_plp; *_tid = iter->tid; *_pos = iter->pos;
		// update iter->tid and iter->pos
		if (iter->head->next) {
			if (iter->tid > iter->head->b.core.tid) {
				fprintf(stderr, "[%s] unsorted input. Pileup aborts.\n", __func__);
				*_n_plp = -1;
				return 0;
			}
		}
		if (iter->tid < iter->head->b.core.tid) { // come to a new reference sequence
			iter->tid = iter->head->b.core.tid; iter->pos = iter->head->beg; // jump to the next reference
		} else if (iter->pos < iter->head->beg) { // here: tid == head->b.core.tid
			iter->pos = iter->head->beg; // jump to the next position
		} else ++iter->pos; // scan contiguously
		// return
		if (n_plp) return iter->plp;
		if (iter->is_eof && iter->head->next == 0) break;
	}
	return 0;
}

int bam_plp_push(bam_plp_t iter, const bam1_t *b)
{
	if (b) {
		if (b->core.tid < 0) return 0;
		if (b->core.flag & iter->flag_mask) return 0;
		bam_copy1(&iter->tail->b, b);
		iter->tail->beg = b->core.pos; iter->tail->end = bam_calend(&b->core, bam1_cigar(b));
		if (b->core.tid < iter->max_tid) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (chromosomes out of order)\n");
			return -1;
		}
		if ((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos)) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (reads out of order)\n");
			return -1;
		}
		iter->max_tid = b->core.tid; iter->max_pos = iter->tail->beg;
		if (iter->tail->end > iter->pos || iter->tail->b.core.tid > iter->tid) {
			iter->tail->next = mp_alloc(iter->mp);
			iter->tail = iter->tail->next;
		}
	} else iter->is_eof = 1;
	return 0;
}

void bam_plp_reset(bam_plp_t iter)
{
	lbnode_t *p, *q;
	iter->max_tid = iter->max_pos = -1;
	iter->tid = iter->pos = 0;
	iter->is_eof = 0;
	for (p = iter->head; p->next;) {
		q = p->next;
		mp_free(iter->mp, p);
		p = q;
	}
	iter->head = iter->tail;
}

void bam_plp_set_mask(bam_plp_t iter, int mask)
{
	iter->flag_mask = mask < 0? BAM_DEF_MASK : (BAM_FUNMAP | mask);
}

void bam_plp_destroy(bam_plp_t iter)
{
	mp_free(iter->mp, iter->dummy);
	mp_free(iter->mp, iter->head);
	if (iter->mp->cnt != 0)
		fprintf(stderr, "[bam_plp_destroy] memory leak: %d. Continue anyway.\n", iter->mp->cnt);
	mp_destroy(iter->mp);
	free(iter->plp);
	free(iter);
}

/*****************
 * callback APIs *
 *****************/

int bam_pileup_file(bamFile fp, int mask, bam_pileup_f func, void *func_data)
{
	bam_plbuf_t *buf;
	int ret;
	bam1_t *b;
	b = bam_init1();
	buf = bam_plbuf_init(func, func_data);
	bam_plbuf_set_mask(buf, mask);
	while ((ret = bam_read1(fp, b)) >= 0)
		bam_plbuf_push(b, buf);
	bam_plbuf_push(0, buf);
	bam_plbuf_destroy(buf);
	bam_destroy1(b);
	return 0;
}

void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask)
{
	bam_plp_set_mask(buf->iter, mask);
}

void bam_plbuf_reset(bam_plbuf_t *buf)
{
	bam_plp_reset(buf->iter);
}

bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)
{
	bam_plbuf_t *buf;
	buf = calloc(1, sizeof(bam_plbuf_t));
	buf->iter = bam_plp_init();
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
	while ((plp = bam_plp_next(buf->iter, &n_plp, &tid, &pos)) != 0)
		buf->func(tid, pos, n_plp, plp, buf->data);
	return 0;
}
