#include <stdio.h>
#include <ctype.h>
#include "bam.h"
#include "bam_endian.h"

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
	assert(header->n_targets > 0);
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

int bam_header_write(bamFile fp, const bam_header_t *header)
{
	char buf[4];
	int32_t i, name_len, x;
	// write "BAM1"
	strncpy(buf, "BAM\001", 4);
	bam_write(fp, buf, 4);
	// write plain text and the number of reference sequences
	if (bam_is_be) {
		x = bam_swap_endian_4(header->l_text);
		bam_write(fp, &x, 4);
		if (header->l_text) bam_write(fp, header->text, header->l_text);
		x = bam_swap_endian_4(header->n_targets);
		bam_write(fp, &x, 4);
	} else {
		bam_write(fp, &header->l_text, 4);
		if (header->l_text) bam_write(fp, header->text, header->l_text);
		bam_write(fp, &header->n_targets, 4);
	}
	// write sequence names and lengths
	for (i = 0; i != header->n_targets; ++i) {
		char *p = header->target_name[i];
		name_len = strlen(p) + 1;
		if (bam_is_be) {
			x = bam_swap_endian_4(name_len);
			bam_write(fp, &x, 4);
		} else bam_write(fp, &name_len, 4);
		bam_write(fp, p, name_len);
		if (bam_is_be) {
			x = bam_swap_endian_4(header->target_len[i]);
			bam_write(fp, &x, 4);
		} else bam_write(fp, &header->target_len[i], 4);
	}
	return 0;
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

inline int bam_write1_core(bamFile fp, const bam1_core_t *c, int data_len, uint8_t *data)
{
	uint32_t x[8], block_len = data_len + BAM_CORE_SIZE, y;
	int i;
	assert(BAM_CORE_SIZE == 32);
	x[0] = c->tid;
	x[1] = c->pos;
	x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | c->l_qname;
	x[3] = (uint32_t)c->flag<<16 | c->n_cigar;
	x[4] = c->l_qseq;
	x[5] = c->mtid;
	x[6] = c->mpos;
	x[7] = c->isize;
	if (bam_is_be) {
		for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
		y = block_len;
		bam_write(fp, bam_swap_endian_4p(&y), 4);
		swap_endian_data(c, data_len, data);
	} else bam_write(fp, &block_len, 4);
	bam_write(fp, x, BAM_CORE_SIZE);
	bam_write(fp, data, data_len);
	if (bam_is_be) swap_endian_data(c, data_len, data);
	return 4 + block_len;
}

int bam_write1(bamFile fp, const bam1_t *b)
{
	return bam_write1_core(fp, &b->core, b->data_len, b->data);
}

void bam_view1(const bam_header_t *header, const bam1_t *b)
{
	uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
	int i;
	const bam1_core_t *c = &b->core;
	printf("%s\t%d\t", bam1_qname(b), c->flag);
	if (c->tid < 0) printf("*\t");
	else printf("%s\t", header->target_name[c->tid]);
	printf("%d\t%d\t", c->pos + 1, c->qual);
	if (c->n_cigar == 0) putchar('*');
	else {
		for (i = 0; i < c->n_cigar; ++i)
			printf("%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	}
	putchar('\t');
	if (c->mtid < 0) printf("*\t");
	else if (c->mtid == c->tid) printf("=\t");
	else printf("%s\t", header->target_name[c->mtid]);
	printf("%d\t%d\t", c->mpos + 1, c->isize);
	for (i = 0; i < c->l_qseq; ++i) putchar(bam_nt16_rev_table[bam1_seqi(s, i)]);
	putchar('\t');
	for (i = 0; i < c->l_qseq; ++i) putchar(t[i] + 33);
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s; ++s;
		printf("\t%c%c:", key[0], key[1]);
		if (type == 'A') { printf("A:%c", *s); ++s; }
		else if (type == 'C') { printf("i:%u", *s); ++s; }
		else if (type == 'c') { printf("i:%d", *s); ++s; }
		else if (type == 'S') { printf("i:%u", *(uint16_t*)s); s += 2; }
		else if (type == 's') { printf("i:%d", *(int16_t*)s); s += 2; }
		else if (type == 'I') { printf("i:%u", *(uint32_t*)s); s += 4; }
		else if (type == 'i') { printf("i:%d", *(int32_t*)s); s += 4; }
		else if (type == 'f') { printf("f:%g", *(float*)s); s += 4; }
		else if (type == 'Z' || type == 'H') { printf("%c:", type); while (*s) putchar(*s++); ++s; }
	}
	putchar('\n');
}
