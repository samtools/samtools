#include <string.h>
#include <stdlib.h>
#include "glf.h"

#ifdef _NO_BGZF
// then alias bgzf_*() functions
#endif

static int glf3_is_BE = 0;

static inline uint32_t bam_swap_endian_4(uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}

static inline uint16_t bam_swap_endian_2(uint16_t v)
{
	return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}

static inline int bam_is_big_endian()
{
	long one= 1;
	return !(*((char *)(&one)));
}

glf3_header_t *glf3_header_init()
{
	glf3_is_BE = bam_is_big_endian();
	return (glf3_header_t*)calloc(1, sizeof(glf3_header_t));
}

glf3_header_t *glf3_header_read(glfFile fp)
{
	glf3_header_t *h;
	char magic[4];
	h = glf3_header_init();
	bgzf_read(fp, magic, 4);
	if (strncmp(magic, "GLF\3", 4)) {
		fprintf(stderr, "[glf3_header_read] invalid magic.\n");
		glf3_header_destroy(h);
		return 0;
	}
	bgzf_read(fp, &h->l_text, 4);
	if (glf3_is_BE) h->l_text = bam_swap_endian_4(h->l_text);
	if (h->l_text) {
		h->text = (uint8_t*)calloc(h->l_text + 1, 1);
		bgzf_read(fp, h->text, h->l_text);
	}
	return h;
}

void glf3_header_write(glfFile fp, const glf3_header_t *h)
{
	int32_t x;
	bgzf_write(fp, "GLF\3", 4);
	x = glf3_is_BE? bam_swap_endian_4(h->l_text) : h->l_text;
	bgzf_write(fp, &x, 4);
	if (h->l_text) bgzf_write(fp, h->text, h->l_text);
}

void glf3_header_destroy(glf3_header_t *h)
{
	free(h->text);
	free(h);
}

char *glf3_ref_read(glfFile fp, int *len)
{
	int32_t n, x;
	char *str;
	*len = 0;
	if (bgzf_read(fp, &n, 4) != 4) return 0;
	if (glf3_is_BE) n = bam_swap_endian_4(n);
	if (n < 0) {
		fprintf(stderr, "[glf3_ref_read] invalid reference name length: %d.\n", n);
		return 0;
	}
	str = (char*)calloc(n + 1, 1); // not necesarily n+1 in fact
	x = bgzf_read(fp, str, n);
	x += bgzf_read(fp, len, 4);
	if (x != n + 4) {
		free(str); *len = -1; return 0; // truncated
	}
	if (glf3_is_BE) *len = bam_swap_endian_4(*len);
	return str;
}

void glf3_ref_write(glfFile fp, const char *str, int len)
{
	int32_t m, n = strlen(str) + 1;
	m = glf3_is_BE? bam_swap_endian_4(n) : n;
	bgzf_write(fp, &m, 4);
	bgzf_write(fp, str, n);
	if (glf3_is_BE) len = bam_swap_endian_4(len);
	bgzf_write(fp, &len, 4);
}

void glf3_view1(const char *ref_name, const glf3_t *g3, int pos)
{
	int j;
	if (g3->rtype == GLF3_RTYPE_END) return;
	printf("%s\t%d\t%c\t%d\t%d\t%d", ref_name, pos + 1,
		   g3->rtype == GLF3_RTYPE_INDEL? '*' : "XACMGRSVTWYHKDBN"[g3->ref_base],
		   g3->depth, g3->rms_mapQ, g3->min_lk);
	if (g3->rtype == GLF3_RTYPE_SUB)
		for (j = 0; j != 10; ++j) printf("\t%d", g3->lk[j]);
	else {
		printf("\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t", g3->lk[0], g3->lk[1], g3->lk[2], g3->indel_len[0], g3->indel_len[1],
			   g3->indel_len[0]? g3->indel_seq[0] : "*", g3->indel_len[1]? g3->indel_seq[1] : "*");
	}
	printf("\n");
}

int glf3_write1(glfFile fp, const glf3_t *g3)
{
	int r;
	uint8_t c;
	uint32_t y[2];
	c = g3->rtype<<4 | g3->ref_base;
	r = bgzf_write(fp, &c, 1);
	if (g3->rtype == GLF3_RTYPE_END) return r;
	y[0] = g3->offset;
	y[1] = g3->min_lk<<24 | g3->depth;
	if (glf3_is_BE) {
		y[0] = bam_swap_endian_4(y[0]);
		y[1] = bam_swap_endian_4(y[1]);
	}
	r += bgzf_write(fp, y, 8);
	r += bgzf_write(fp, &g3->rms_mapQ, 1);
	if (g3->rtype == GLF3_RTYPE_SUB) r += bgzf_write(fp, g3->lk, 10);
	else {
		int16_t x[2];
		r += bgzf_write(fp, g3->lk, 3);
		x[0] = glf3_is_BE? bam_swap_endian_2(g3->indel_len[0]) : g3->indel_len[0];
		x[1] = glf3_is_BE? bam_swap_endian_2(g3->indel_len[1]) : g3->indel_len[1];
		r += bgzf_write(fp, x, 4);
		if (g3->indel_len[0]) r += bgzf_write(fp, g3->indel_seq[0], abs(g3->indel_len[0]));
		if (g3->indel_len[1]) r += bgzf_write(fp, g3->indel_seq[1], abs(g3->indel_len[1]));
	}
	return r;
}

#ifndef kv_roundup32
#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

int glf3_read1(glfFile fp, glf3_t *g3)
{
	int r;
	uint8_t c;
	uint32_t y[2];
	r = bgzf_read(fp, &c, 1);
	if (r == 0) return 0;
	g3->ref_base = c & 0xf;
	g3->rtype = c>>4;
	if (g3->rtype == GLF3_RTYPE_END) return r;
	r += bgzf_read(fp, y, 8);
	if (glf3_is_BE) {
		y[0] = bam_swap_endian_4(y[0]);
		y[1] = bam_swap_endian_4(y[1]);
	}
	g3->offset = y[0];
	g3->min_lk = y[1]>>24;
	g3->depth = y[1]<<8>>8;
	r += bgzf_read(fp, &g3->rms_mapQ, 1);
	if (g3->rtype == GLF3_RTYPE_SUB) r += bgzf_read(fp, g3->lk, 10);
	else {
		int16_t x[2], max;
		r += bgzf_read(fp, g3->lk, 3);
		r += bgzf_read(fp, x, 4);
		if (glf3_is_BE) {
			x[0] = bam_swap_endian_2(x[0]);
			x[1] = bam_swap_endian_2(x[1]);
		}
		g3->indel_len[0] = x[0];
		g3->indel_len[1] = x[1];
		x[0] = abs(x[0]); x[1] = abs(x[1]);
		max = (x[0] > x[1]? x[0] : x[1]) + 1;
		if (g3->max_len < max) {
			g3->max_len = max;
			kv_roundup32(g3->max_len);
			g3->indel_seq[0] = (char*)realloc(g3->indel_seq[0], g3->max_len);
			g3->indel_seq[1] = (char*)realloc(g3->indel_seq[1], g3->max_len);
		}
		r += bgzf_read(fp, g3->indel_seq[0], x[0]);
		r += bgzf_read(fp, g3->indel_seq[1], x[1]);
		g3->indel_seq[0][x[0]] = g3->indel_seq[1][x[1]] = 0;
	}
	return r;
}

void glf3_view(glfFile fp)
{
	glf3_header_t *h;
	char *name;
	glf3_t *g3;
	int len;
	h = glf3_header_read(fp);
	g3 = glf3_init1();
	while ((name = glf3_ref_read(fp, &len)) != 0) {
		int pos = 0;
		while (glf3_read1(fp, g3) && g3->rtype != GLF3_RTYPE_END) {
			pos += g3->offset;
			glf3_view1(name, g3, pos);
		}
		free(name);
	}
	glf3_header_destroy(h);
	glf3_destroy1(g3);
}

int glf3_view_main(int argc, char *argv[])
{
	glfFile fp;
	if (argc == 1) {
		fprintf(stderr, "Usage: glfview <in.glf>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? bgzf_fdopen(fileno(stdin), "r") : bgzf_open(argv[1], "r");
	if (fp == 0) {
		fprintf(stderr, "Fail to open file '%s'\n", argv[1]);
		return 1;
	}
	glf3_view(fp);
	bgzf_close(fp);
	return 0;
}

#ifdef GLFVIEW_MAIN
int main(int argc, char *argv[])
{
	return glf3_view_main(argc, argv);
}
#endif
