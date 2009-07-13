#include <ctype.h>
#include <assert.h>
#include <stdarg.h>

#define BAM_LITE
#include "bam.h"

#undef bam_open
#undef bam_dopen
#undef bam_close
#undef bam_read
#undef bam_write
#undef bam_tell
#undef bam_seek

#define bam_open(fn, mode) gzopen(fn, mode)
#define bam_dopen(fd, mode) gzdopen(fd, mode)
#define bam_close(fp) gzclose(fp)
#define bam_read(fp, buf, size) gzread(fp, buf, size)

#define BAM_NO_HASH

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

/*************************
 *** from bam_endian.h ***
 *************************/

static inline int bam_is_big_endian()
{
	long one= 1;
	return !(*((char *)(&one)));
}
static inline uint16_t bam_swap_endian_2(uint16_t v)
{
	return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *bam_swap_endian_2p(void *x)
{
	*(uint16_t*)x = bam_swap_endian_2(*(uint16_t*)x);
	return x;
}
static inline uint32_t bam_swap_endian_4(uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *bam_swap_endian_4p(void *x)
{
	*(uint32_t*)x = bam_swap_endian_4(*(uint32_t*)x);
	return x;
}
static inline uint64_t bam_swap_endian_8(uint64_t v)
{
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *bam_swap_endian_8p(void *x)
{
	*(uint64_t*)x = bam_swap_endian_8(*(uint64_t*)x);
	return x;
}

/**********************
 *** from kstring.* ***
 **********************/

typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	strncpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

static inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

static inline int kputc(int c, kstring_t *s)
{
	if (s->l + 1 >= s->m) {
		s->m = s->l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	s->s[s->l++] = c;
	s->s[s->l] = 0;
	return c;
}

int ksprintf(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	int l;
	va_start(ap, fmt);
	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap); // This line does not work with glibc 2.0. See `man snprintf'.
	va_end(ap);
	if (l + 1 > s->m - s->l) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
		va_start(ap, fmt);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap);
	}
	va_end(ap);
	s->l += l;
	return l;
}

/******************
 *** from bam.c ***
 ******************/

int bam_is_be = 0;

/**************************
 * CIGAR related routines *
 **************************/

int bam_segreg(int32_t pos, const bam1_core_t *c, const uint32_t *cigar, bam_segreg_t *reg)
{
	unsigned k;
	int32_t x = c->pos, y = 0;
	int state = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK; // operation
		int l = cigar[k] >> BAM_CIGAR_SHIFT; // length
		if (state == 0 && (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CINS) && x + l > pos) {
			reg->tbeg = x; reg->qbeg = y; reg->cbeg = k;
			state = 1;
		}
		if (op == BAM_CMATCH) { x += l; y += l; }
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP) x += l;
		else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		if (state == 1 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CREF_SKIP || k == c->n_cigar - 1)) {
			reg->tend = x; reg->qend = y; reg->cend = k;
		}
	}
	return state? 0 : -1;
}

uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar)
{
	uint32_t k, end;
	end = c->pos;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
			end += cigar[k] >> BAM_CIGAR_SHIFT;
	}
	return end;
}

int32_t bam_cigar2qlen(const bam1_core_t *c, const uint32_t *cigar)
{
	uint32_t k;
	int32_t l = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP)
			l += cigar[k] >> BAM_CIGAR_SHIFT;
	}
	return l;
}

/********************
 * BAM I/O routines *
 ********************/

bam_header_t *bam_header_init()
{
	bam_is_be = bam_is_big_endian();
	return (bam_header_t*)calloc(1, sizeof(bam_header_t));
}

void bam_header_destroy(bam_header_t *header)
{
	int32_t i;
	extern void bam_destroy_header_hash(bam_header_t *header);
	if (header == 0) return;
	if (header->target_name) {
		for (i = 0; i < header->n_targets; ++i)
			free(header->target_name[i]);
		free(header->target_name);
		free(header->target_len);
	}
	free(header->text);
#ifndef BAM_NO_HASH
	if (header->rg2lib) bam_strmap_destroy(header->rg2lib);
	bam_destroy_header_hash(header);
#endif
	free(header);
}

bam_header_t *bam_header_read(bamFile fp)
{
	bam_header_t *header;
	char buf[4];
	int32_t i, name_len;
	// read "BAM1"
	if (bam_read(fp, buf, 4) != 4) return 0;
	if (strncmp(buf, "BAM\001", 4)) {
		fprintf(stderr, "[bam_header_read] wrong header\n");
		return 0;
	}
	header = bam_header_init();
	// read plain text and the number of reference sequences
	bam_read(fp, &header->l_text, 4);
	if (bam_is_be) bam_swap_endian_4p(&header->l_text);
	header->text = (char*)calloc(header->l_text + 1, 1);
	bam_read(fp, header->text, header->l_text);
	bam_read(fp, &header->n_targets, 4);
	if (bam_is_be) bam_swap_endian_4p(&header->n_targets);
	// read reference sequence names and lengths
	header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
	header->target_len = (uint32_t*)calloc(header->n_targets, 4);
	for (i = 0; i != header->n_targets; ++i) {
		bam_read(fp, &name_len, 4);
		if (bam_is_be) bam_swap_endian_4p(&name_len);
		header->target_name[i] = (char*)calloc(name_len, 1);
		bam_read(fp, header->target_name[i], name_len);
		bam_read(fp, &header->target_len[i], 4);
		if (bam_is_be) bam_swap_endian_4p(&header->target_len[i]);
	}
	return header;
}

static void swap_endian_data(const bam1_core_t *c, int data_len, uint8_t *data)
{
	uint8_t *s;
	uint32_t i, *cigar = (uint32_t*)(data + c->l_qname);
	s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
	for (i = 0; i < c->n_cigar; ++i) bam_swap_endian_4p(&cigar[i]);
	while (s < data + data_len) {
		uint8_t type;
		s += 2; // skip key
		type = toupper(*s); ++s; // skip type
		if (type == 'C' || type == 'A') ++s;
		else if (type == 'S') { bam_swap_endian_2p(s); s += 2; }
		else if (type == 'I' || type == 'F') { bam_swap_endian_4p(s); s += 4; }
		else if (type == 'D') { bam_swap_endian_8p(s); s += 8; }
		else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
	}
}

int bam_read1(bamFile fp, bam1_t *b)
{
	bam1_core_t *c = &b->core;
	int32_t block_len, ret, i;
	uint32_t x[8];

	assert(BAM_CORE_SIZE == 32);
	if ((ret = bam_read(fp, &block_len, 4)) != 4) {
		if (ret == 0) return -1; // normal end-of-file
		else return -2; // truncated
	}
	if (bam_read(fp, x, BAM_CORE_SIZE) != BAM_CORE_SIZE) return -3;
	if (bam_is_be) {
		bam_swap_endian_4p(&block_len);
		for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
	}
	c->tid = x[0]; c->pos = x[1];
	c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
	c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
	b->data_len = block_len - BAM_CORE_SIZE;
	if (b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	if (bam_read(fp, b->data, b->data_len) != b->data_len) return -4;
	b->l_aux = b->data_len - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
	if (bam_is_be) swap_endian_data(c, b->data_len, b->data);
	return 4 + block_len;
}

char *bam_format1(const bam_header_t *header, const bam1_t *b)
{
	uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
	int i;
	const bam1_core_t *c = &b->core;
	kstring_t str;
	str.l = str.m = 0; str.s = 0;

	ksprintf(&str, "%s\t%d\t", bam1_qname(b), c->flag);
	if (c->tid < 0) kputs("*\t", &str);
	else ksprintf(&str, "%s\t", header->target_name[c->tid]);
	ksprintf(&str, "%d\t%d\t", c->pos + 1, c->qual);
	if (c->n_cigar == 0) kputc('*', &str);
	else {
		for (i = 0; i < c->n_cigar; ++i)
			ksprintf(&str, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	}
	kputc('\t', &str);
	if (c->mtid < 0) kputs("*\t", &str);
	else if (c->mtid == c->tid) kputs("=\t", &str);
	else ksprintf(&str, "%s\t", header->target_name[c->mtid]);
	ksprintf(&str, "%d\t%d\t", c->mpos + 1, c->isize);
	for (i = 0; i < c->l_qseq; ++i) kputc(bam_nt16_rev_table[bam1_seqi(s, i)], &str);
	kputc('\t', &str);
	if (t[0] == 0xff) kputc('*', &str);
	else for (i = 0; i < c->l_qseq; ++i) kputc(t[i] + 33, &str);
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s; ++s;
		ksprintf(&str, "\t%c%c:", key[0], key[1]);
		if (type == 'A') { ksprintf(&str, "A:%c", *s); ++s; }
		else if (type == 'C') { ksprintf(&str, "i:%u", *s); ++s; }
		else if (type == 'c') { ksprintf(&str, "i:%d", *s); ++s; }
		else if (type == 'S') { ksprintf(&str, "i:%u", *(uint16_t*)s); s += 2; }
		else if (type == 's') { ksprintf(&str, "i:%d", *(int16_t*)s); s += 2; }
		else if (type == 'I') { ksprintf(&str, "i:%u", *(uint32_t*)s); s += 4; }
		else if (type == 'i') { ksprintf(&str, "i:%d", *(int32_t*)s); s += 4; }
		else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
		else if (type == 'd') { ksprintf(&str, "d:%lg", *(double*)s); s += 8; }
		else if (type == 'Z' || type == 'H') { ksprintf(&str, "%c:", type); while (*s) kputc(*s++, &str); ++s; }
	}
	return str.s;
}

/*************************
 *** from bam_pileup.c ***
 *************************/

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

struct __bam_plbuf_t {
	mempool_t *mp;
	lbnode_t *head, *tail, *dummy;
	bam_pileup_f func;
	void *func_data;
	int32_t tid, pos, max_tid, max_pos;
	int max_pu, is_eof;
	bam_pileup1_t *pu;
	int flag_mask;
};

void bam_plbuf_reset(bam_plbuf_t *buf)
{
	lbnode_t *p, *q;
	buf->max_tid = buf->max_pos = -1;
	buf->tid = buf->pos = 0;
	buf->is_eof = 0;
	for (p = buf->head; p->next;) {
		q = p->next;
		mp_free(buf->mp, p);
		p = q;
	}
	buf->head = buf->tail;
}

void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask)
{
	if (mask < 0) buf->flag_mask = BAM_DEF_MASK;
	else buf->flag_mask = BAM_FUNMAP | mask;
}

bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)
{
	bam_plbuf_t *buf;
	buf = (bam_plbuf_t*)calloc(1, sizeof(bam_plbuf_t));
	buf->func = func; buf->func_data = data;
	buf->mp = mp_init();
	buf->head = buf->tail = mp_alloc(buf->mp);
	buf->dummy = mp_alloc(buf->mp);
	buf->max_tid = buf->max_pos = -1;
	buf->flag_mask = BAM_DEF_MASK;
	return buf;
}

void bam_plbuf_destroy(bam_plbuf_t *buf)
{
	mp_free(buf->mp, buf->dummy);
	mp_free(buf->mp, buf->head);
	if (buf->mp->cnt != 0)
		fprintf(stderr, "[bam_plbuf_destroy] memory leak: %d. Continue anyway.\n", buf->mp->cnt);
	mp_destroy(buf->mp);
	free(buf->pu);
	free(buf);
}

int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf)
{
	if (b) { // fill buffer
		if (b->core.tid < 0) return 0;
		if (b->core.flag & buf->flag_mask) return 0;
		bam_copy1(&buf->tail->b, b);
		buf->tail->beg = b->core.pos; buf->tail->end = bam_calend(&b->core, bam1_cigar(b));
		if (!(b->core.tid >= buf->max_tid || (b->core.tid == buf->max_tid && buf->tail->beg >= buf->max_pos))) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted. Abort!\n");
			abort();
		}
		buf->max_tid = b->core.tid; buf->max_pos = buf->tail->beg;
		if (buf->tail->end > buf->pos || buf->tail->b.core.tid > buf->tid) {
			buf->tail->next = mp_alloc(buf->mp);
			buf->tail = buf->tail->next;
		}
	} else buf->is_eof = 1;
	while (buf->is_eof || buf->max_tid > buf->tid || (buf->max_tid == buf->tid && buf->max_pos > buf->pos)) {
		int n_pu = 0;
		lbnode_t *p, *q;
		buf->dummy->next = buf->head;
		for (p = buf->head, q = buf->dummy; p->next; q = p, p = p->next) {
			if (p->b.core.tid < buf->tid || (p->b.core.tid == buf->tid && p->end <= buf->pos)) { // then remove from the list
				q->next = p->next; mp_free(buf->mp, p); p = q;
			} else if (p->b.core.tid == buf->tid && p->beg <= buf->pos) { // here: p->end > pos; then add to pileup
				if (n_pu == buf->max_pu) { // then double the capacity
					buf->max_pu = buf->max_pu? buf->max_pu<<1 : 256;
					buf->pu = (bam_pileup1_t*)realloc(buf->pu, sizeof(bam_pileup1_t) * buf->max_pu);
				}
				buf->pu[n_pu].b = &p->b;
				if (resolve_cigar(buf->pu + n_pu, buf->pos)) ++n_pu; // skip the read if we are looking at BAM_CREF_SKIP
			}
		}
		buf->head = buf->dummy->next; // dummy->next may be changed
		if (n_pu) { // then call user defined function
			buf->func(buf->tid, buf->pos, n_pu, buf->pu, buf->func_data);
		}
		// update tid and pos
		if (buf->head->next) {
			if (buf->tid > buf->head->b.core.tid) {
				fprintf(stderr, "[bam_plbuf_push] unsorted input. Pileup aborts.\n");
				return 1;
			}
		}
		if (buf->tid < buf->head->b.core.tid) { // come to a new reference sequence
			buf->tid = buf->head->b.core.tid; buf->pos = buf->head->beg; // jump to the next reference
		} else if (buf->pos < buf->head->beg) { // here: tid == head->b.core.tid
			buf->pos = buf->head->beg; // jump to the next position
		} else ++buf->pos; // scan contiguously
		if (buf->is_eof && buf->head->next == 0) break;
	}
	return 0;
}

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

/**********************
 *** from bam_aux.c ***
 **********************/

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		int type, x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return s;
		type = toupper(*s); ++s;
		if (type == 'C') ++s;
		else if (type == 'S') s += 2;
		else if (type == 'I' || type == 'F') s += 4;
		else if (type == 'D') s += 8;
		else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
	}
	return 0;
}
