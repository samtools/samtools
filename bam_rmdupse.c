#include <math.h>
#include "sam.h"
#include "khash.h"

typedef struct {
	int n, m;
	int *a;
} listelem_t;

KHASH_MAP_INIT_INT(32, listelem_t)

//#define BLOCK_SIZE 65536
#define BLOCK_SIZE 100

typedef struct {
	bam1_t *b;
	int rpos, score;
} elem_t;

typedef struct {
	int n, max, x;
	elem_t *buf;
} buffer_t;

static int fill_buf(samfile_t *in, buffer_t *buf)
{
	int i, ret, last_tid, min_rpos = 0x7fffffff, capacity;
	bam1_t *b = bam_init1();
	bam1_core_t *c = &b->core;
	// squeeze out the empty cells at the beginning
	for (i = 0; i < buf->n; ++i)
		if (buf->buf[i].b) break;
	if (i < buf->n) { // squeeze
		if (i > 0) {
			memmove(buf->buf, buf->buf + i, sizeof(elem_t) * (buf->n - i));
			buf->n = buf->n - i;
		}
	} else buf->n = 0;
	// calculate min_rpos
	for (i = 0; i < buf->n; ++i) {
		elem_t *e = buf->buf + i;
		if (e->b && e->rpos >= 0 && e->rpos < min_rpos)
			min_rpos = buf->buf[i].rpos;
	}
	// fill the buffer
	buf->x = -1;
	last_tid = buf->n? buf->buf[0].b->core.tid : -1;
	capacity = buf->n + BLOCK_SIZE;
	while ((ret = samread(in, b)) >= 0) {
		elem_t *e;
		uint8_t *qual = bam1_qual(b);
		int is_mapped;
		if (last_tid < 0) last_tid = c->tid;
		if (c->tid != last_tid) {
			if (buf->x < 0) buf->x = buf->n;
		}
		if (buf->n >= buf->max) { // enlarge
			buf->max = buf->max? buf->max<<1 : 8;
			buf->buf = (elem_t*)realloc(buf->buf, sizeof(elem_t) * buf->max);
		}
		e = &buf->buf[buf->n++];
		e->b = bam_dup1(b);
		e->rpos = -1; e->score = 0;
		for (i = 0; i < c->l_qseq; ++i) e->score += qual[i] + 1;
		e->score = (double)e->score / sqrt(c->l_qseq + 1);
		is_mapped = (c->tid < 0 || c->tid >= in->header->n_targets || (c->flag&BAM_FUNMAP))? 0 : 1;
		if (is_mapped && (c->flag & BAM_FREVERSE)) {
			e->rpos = b->core.pos + bam_calend(&b->core, bam1_cigar(b));
			if (min_rpos > e->rpos) min_rpos = e->rpos;
		}
		if (buf->n >= capacity) {
			if (c->pos <= min_rpos) capacity += BLOCK_SIZE;
			else break;
		}
	}
	if (ret >= 0 && buf->x < 0) buf->x = buf->n;
	bam_destroy1(b);
	return buf->n;
}

static void rmdupse_buf(buffer_t *buf)
{
	khash_t(32) *h;
	uint32_t key;
	khint_t k;
	int mpos, i, upper;
	listelem_t *p;
	mpos = 0x7fffffff;
	mpos = (buf->x == buf->n)? buf->buf[buf->x-1].b->core.pos : 0x7fffffff;
	upper = (buf->x < 0)? buf->n : buf->x;
	// fill the hash table
	h = kh_init(32);
	for (i = 0; i < upper; ++i) {
		elem_t *e = buf->buf + i;
		int ret;
		if (e->score < 0) continue;
		if (e->rpos >= 0) {
			if (e->rpos <= mpos) key = (uint32_t)e->rpos<<1 | 1;
			else continue;
		} else {
			if (e->b->core.pos < mpos) key = (uint32_t)e->b->core.pos<<1;
			else continue;
		}
		k = kh_put(32, h, key, &ret);
		p = &kh_val(h, k);
		if (ret == 0) { // present in the hash table
			if (p->n == p->m) {
				p->m <<= 1;
				p->a = (int*)realloc(p->a, p->m * sizeof(int));
			}
			p->a[p->n++] = i;
		} else {
			p->m = p->n = 1;
			p->a = (int*)calloc(p->m, sizeof(int));
			p->a[0] = i;
		}
	}
	// rmdup
	for (k = kh_begin(h); k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			int max, maxi;
			p = &kh_val(h, k);
			// get the max
			for (i = max = 0, maxi = -1; i < p->n; ++i) {
				if (buf->buf[p->a[i]].score > max) {
					max = buf->buf[p->a[i]].score;
					maxi = i;
				}
			}
			// mark the elements
			for (i = 0; i < p->n; ++i) {
				buf->buf[p->a[i]].score = -1;
				if (i != maxi) {
					bam_destroy1(buf->buf[p->a[i]].b);
					buf->buf[p->a[i]].b = 0;
				}
			}
			// free
			free(p->a);
		}
	}
	kh_destroy(32, h);
}

static void dump_buf(buffer_t *buf, samfile_t *out)
{
	int i;
	for (i = 0; i < buf->n; ++i) {
		elem_t *e = buf->buf + i;
		if (e->score != -1) break;
		if (e->b) {
			samwrite(out, e->b);
			bam_destroy1(e->b);
			e->b = 0;
		}
	}
}

int bam_rmdupse(int argc, char *argv[])
{
	samfile_t *in, *out;
	buffer_t *buf;
	if (argc < 3) {
		fprintf(stderr, "Usage: samtools rmdupse <in.bam> <out.bam>\n");
		return 1;
	}
	buf = calloc(1, sizeof(buffer_t));
	in = samopen(argv[1], "rb", 0);
	out = samopen(argv[2], "wb", in->header);
	while (fill_buf(in, buf)) {
		rmdupse_buf(buf);
		dump_buf(buf, out);
	}
	samclose(in); samclose(out);
	free(buf->buf); free(buf);
	return 0;
}
