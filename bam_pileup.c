#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "sam.h"

typedef struct {
	int k, x, y, end;
} cstate_t;

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

typedef struct __linkbuf_t {
	bam1_t b;
	uint32_t beg, end;
	cstate_t s;
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

/* s->k: the index of the CIGAR operator that has just been processed.
   s->x: the reference coordinate of the start of s->k
   s->y: the query coordiante of the start of s->k
 */
static inline int resolve_cigar2(bam_pileup1_t *p, uint32_t pos, cstate_t *s)
{
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)

	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t *cigar = bam1_cigar(b);
	int k, is_head = 0;
	// determine the current CIGAR operation
//	fprintf(stderr, "%s\tpos=%d\tend=%d\t(%d,%d,%d)\n", bam1_qname(b), pos, s->end, s->k, s->x, s->y);
	if (s->k == -1) { // never processed
		is_head = 1;
		if (c->n_cigar == 1) { // just one operation, save a loop
			if (_cop(cigar[0]) == BAM_CMATCH) s->k = 0, s->x = c->pos, s->y = 0;
		} else { // find the first match or deletion
			for (k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
				int op = _cop(cigar[k]);
				int l = _cln(cigar[k]);
				if (op == BAM_CMATCH || op == BAM_CDEL) break;
				else if (op == BAM_CREF_SKIP) s->x += l;
				else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
			}
			assert(k < c->n_cigar);
			s->k = k;
		}
	} else { // the read has been processed before
		int op, l = _cln(cigar[s->k]);
		if (pos - s->x >= l) { // jump to the next operation
			assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
			op = _cop(cigar[s->k+1]);
			if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP) { // jump to the next without a loop
				if (_cop(cigar[s->k]) == BAM_CMATCH) s->y += l;
				s->x += l;
				++s->k;
			} else { // find the next M/D/N
				if (_cop(cigar[s->k]) == BAM_CMATCH) s->y += l;
				s->x += l;
				for (k = s->k + 1; k < c->n_cigar; ++k) {
					op = _cop(cigar[k]), l = _cln(cigar[k]);
					if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP) break;
					else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
				}
				s->k = k;
			}
			assert(s->k < c->n_cigar); // otherwise a bug
		} // else, do nothing
	}
	{ // collect pileup information
		int op, l;
		op = _cop(cigar[s->k]); l = _cln(cigar[s->k]);
		p->is_del = p->indel = p->is_refskip = 0;
		if (s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
			int op2 = _cop(cigar[s->k+1]);
			int l2 = _cln(cigar[s->k+1]);
			if (op2 == BAM_CDEL) p->indel = -(int)l2;
			else if (op2 == BAM_CINS) p->indel = l2;
			else if (op2 == BAM_CPAD && s->k + 2 < c->n_cigar) { // no working for adjacent padding
				int l3 = 0;
				for (k = s->k + 2; k < c->n_cigar; ++k) {
					op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
					if (op2 == BAM_CINS) l3 += l2;
					else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP) break;
				}
				if (l3 > 0) p->indel = l3;
			}
		}
		if (op == BAM_CMATCH) {
			p->qpos = s->y + (pos - s->x);
		} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
			p->is_del = 1; p->qpos = s->y; // FIXME: distinguish D and N!!!!!
			p->is_refskip = (op == BAM_CREF_SKIP);
		} // cannot be other operations; otherwise a bug
		p->is_head = (pos == c->pos); p->is_tail = (pos == s->end);
	}
	return 1;
}

/* --- END: Auxiliary functions */

/*******************
 * pileup iterator *
 *******************/

struct __bam_plp_t {
	mempool_t *mp;
	lbnode_t *head, *tail, *dummy;
	int32_t tid, pos, max_tid, max_pos;
	int is_eof, flag_mask, max_plp, error, maxcnt;
	bam_pileup1_t *plp;
	// for the "auto" interface only
	bam1_t *b;
	bam_plp_auto_f func;
	void *data;
};

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data)
{
	bam_plp_t iter;
	iter = calloc(1, sizeof(struct __bam_plp_t));
	iter->mp = mp_init();
	iter->head = iter->tail = mp_alloc(iter->mp);
	iter->dummy = mp_alloc(iter->mp);
	iter->max_tid = iter->max_pos = -1;
	iter->flag_mask = BAM_DEF_MASK;
	iter->maxcnt = 8000;
	if (func) {
		iter->func = func;
		iter->data = data;
		iter->b = bam_init1();
	}
	return iter;
}

void bam_plp_destroy(bam_plp_t iter)
{
	mp_free(iter->mp, iter->dummy);
	mp_free(iter->mp, iter->head);
	if (iter->mp->cnt != 0)
		fprintf(stderr, "[bam_plp_destroy] memory leak: %d. Continue anyway.\n", iter->mp->cnt);
	mp_destroy(iter->mp);
	if (iter->b) bam_destroy1(iter->b);
	free(iter->plp);
	free(iter);
}

const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
	if (iter->error) { *_n_plp = -1; return 0; }
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
				if (resolve_cigar2(iter->plp + n_plp, iter->pos, &p->s)) ++n_plp; // actually always true...
			}
		}
		iter->head = iter->dummy->next; // dummy->next may be changed
		*_n_plp = n_plp; *_tid = iter->tid; *_pos = iter->pos;
		// update iter->tid and iter->pos
		if (iter->head->next) {
			if (iter->tid > iter->head->b.core.tid) {
				fprintf(stderr, "[%s] unsorted input. Pileup aborts.\n", __func__);
				iter->error = 1;
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
	if (iter->error) return -1;
	if (b) {
		if (b->core.tid < 0) return 0;
		if (b->core.flag & iter->flag_mask) return 0;
		if (iter->tid == b->core.tid && iter->pos == b->core.pos && iter->mp->cnt > iter->maxcnt) return 0;
		bam_copy1(&iter->tail->b, b);
		iter->tail->beg = b->core.pos; iter->tail->end = bam_calend(&b->core, bam1_cigar(b));
		iter->tail->s = g_cstate_null; iter->tail->s.end = iter->tail->end - 1; // initialize cstate_t
		if (b->core.tid < iter->max_tid) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (chromosomes out of order)\n");
			iter->error = 1;
			return -1;
		}
		if ((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos)) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (reads out of order)\n");
			iter->error = 1;
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

const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
	const bam_pileup1_t *plp;
	if (iter->func == 0 || iter->error) { *_n_plp = -1; return 0; }
	if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
	else { // no pileup line can be obtained; read alignments
		*_n_plp = 0;
		if (iter->is_eof) return 0;
		while (iter->func(iter->data, iter->b) >= 0) {
			if (bam_plp_push(iter, iter->b) < 0) {
				*_n_plp = -1;
				return 0;
			}
			if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
			// otherwise no pileup line can be returned; read the next alignment.
		}
		bam_plp_push(iter, 0);
		if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
		return 0;
	}
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

void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)
{
	iter->maxcnt = maxcnt;
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

/***********
 * mpileup *
 ***********/

struct __bam_mplp_t {
	int n;
	uint64_t min, *pos;
	bam_plp_t *iter;
	int *n_plp;
	const bam_pileup1_t **plp;
};

bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data)
{
	int i;
	bam_mplp_t iter;
	iter = calloc(1, sizeof(struct __bam_mplp_t));
	iter->pos = calloc(n, 8);
	iter->n_plp = calloc(n, sizeof(int));
	iter->plp = calloc(n, sizeof(void*));
	iter->iter = calloc(n, sizeof(void*));
	iter->n = n;
	iter->min = (uint64_t)-1;
	for (i = 0; i < n; ++i) {
		iter->iter[i] = bam_plp_init(func, data[i]);
		iter->pos[i] = iter->min;
	}
	return iter;
}

void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt)
{
	int i;
	for (i = 0; i < iter->n; ++i)
		iter->iter[i]->maxcnt = maxcnt;
}

void bam_mplp_destroy(bam_mplp_t iter)
{
	int i;
	for (i = 0; i < iter->n; ++i) bam_plp_destroy(iter->iter[i]);
	free(iter->iter); free(iter->pos); free(iter->n_plp); free(iter->plp);
	free(iter);
}

int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)
{
	int i, ret = 0;
	uint64_t new_min = (uint64_t)-1;
	for (i = 0; i < iter->n; ++i) {
		if (iter->pos[i] == iter->min) {
			int tid, pos;
			iter->plp[i] = bam_plp_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
			iter->pos[i] = (uint64_t)tid<<32 | pos;
		}
		if (iter->plp[i] && iter->pos[i] < new_min) new_min = iter->pos[i];
	}
	iter->min = new_min;
	if (new_min == (uint64_t)-1) return 0;
	*_tid = new_min>>32; *_pos = (uint32_t)new_min;
	for (i = 0; i < iter->n; ++i) {
		if (iter->pos[i] == iter->min) { // FIXME: valgrind reports "uninitialised value(s) at this line"
			n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
			++ret;
		} else n_plp[i] = 0, plp[i] = 0;
	}
	return ret;
}
