#include <errno.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "bcf.h"

static int8_t nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4, -1, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4, -1, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static int read_I16(bcf1_t *b, int anno[16])
{
	char *p;
	int i;
	if ((p = strstr(b->info, "I16=")) == 0) return -1;
	p += 4;
	for (i = 0; i < 16; ++i) {
		anno[i] = strtol(p, &p, 10);
		if (anno[i] == 0 && (errno == EINVAL || errno == ERANGE)) return -2;
		++p;
	}
	return 0;
}

int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b)
{
	int a[4], k, g[10], l, map[4], k1, j, i, i0, anno[16], dp, mq, d_rest;
	char *s;
	if (b->ref[1] != 0 || b->n_alleles > 4) return -1; // ref is not a single base
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) break;
	if (i == b->n_gi) return -1; // no PL
	if (read_I16(b, anno) != 0) return -1; // no I16; FIXME: can be improved
	d_rest = dp = anno[0] + anno[1] + anno[2] + anno[3];
	if (dp == 0) return -1; // depth is zero
	mq = (int)(sqrt((double)(anno[9] + anno[11]) / dp) + .499);
	i0 = i;
	a[0] = nt4_table[(int)b->ref[0]];
	if (a[0] > 3) return -1; // ref is not A/C/G/T
	a[1] = a[2] = a[3] = -2; // -1 has a special meaning
	if (b->alt[0] == 0) return -1; // no alternate allele
	map[0] = map[1] = map[2] = map[3] = -2;
	map[a[0]] = 0;
	for (k = 0, s = b->alt, k1 = -1; k < 3 && *s; ++k, s += 2) {
		if (s[1] != ',' && s[1] != 0) return -1; // ALT is not single base
		a[k+1] = nt4_table[(int)*s];
		if (a[k+1] >= 0) map[a[k+1]] = k+1;
		else k1 = k+1;
		if (s[1] == 0) break;
	}
	for (k = 0; k < 4; ++k)
		if (map[k] < 0) map[k] = k1;
	for (i = 0; i < h->n_smpl; ++i) {
		int d;
		uint8_t *p = b->gi[i0].data + i * b->gi[i0].len;
		for (j = 0; j < b->gi[i0].len; ++j)
			if (p[j]) break;
		d = (int)((double)d_rest / (h->n_smpl - i) + .499);
		if (d == 0) d = 1;
		if (j == b->gi[i0].len) d = 0;
		d_rest -= d;
		for (k = j = 0; k < 4; ++k) {
			for (l = k; l < 4; ++l) {
				int t, x = map[k], y = map[l];
				if (x > y) t = x, x = y, y = t; // swap
				g[j++] = p[y * (y+1) / 2 + x];
			}
		}
		printf("%s\t%d\t%c", h->ns[b->tid], b->pos+1, *b->ref);
		printf("\t%d\t%d\t0", d, mq);
		for (j = 0; j < 10; ++j)
			printf("\t%d", g[j]);
		printf("\t%s\n", h->sns[i]);
	}
	return 0;
}
