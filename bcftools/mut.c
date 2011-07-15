#include <stdlib.h>
#include <stdint.h>
#include "bcf.h"

#define MAX_GENO 359

int8_t seq_bitcnt[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
char *seq_nt16rev = "XACMGRSVTWYHKDBN";

uint32_t *bcf_trio_prep(int is_x, int is_son)
{
	int i, j, k, n, map[10];
	uint32_t *ret;
	ret = calloc(MAX_GENO, 4);
	for (i = 0, k = 0; i < 4; ++i)
		for (j = i; j < 4; ++j)
			map[k++] = 1<<i|1<<j;
	for (i = 0, n = 1; i < 10; ++i) { // father
		if (is_x && seq_bitcnt[map[i]] != 1) continue;
		if (is_x && is_son) {
			for (j = 0; j < 10; ++j) // mother
				for (k = 0; k < 10; ++k) // child
					if (seq_bitcnt[map[k]] == 1 && (map[j]&map[k]))
						ret[n++] = j<<16 | i<<8 | k;
		} else {
			for (j = 0; j < 10; ++j) // mother
				for (k = 0; k < 10; ++k) // child
					if ((map[i]&map[k]) && (map[j]&map[k]) && ((map[i]|map[j])&map[k]) == map[k])
						ret[n++] = j<<16 | i<<8 | k;
		}
	}
	ret[0] = n - 1;
	return ret;
}


int bcf_trio_call(const uint32_t *prep, const bcf1_t *b, int *llr, int64_t *gt)
{
	int i, j, k;
	const bcf_ginfo_t *PL;
	uint8_t *gl10;
	int map[10];
	if (b->n_smpl != 3) return -1; // not a trio
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) break;
	if (i == b->n_gi) return -1; // no PL
	gl10 = alloca(10 * b->n_smpl);
	if (bcf_gl10(b, gl10) < 0) {
		if (bcf_gl10_indel(b, gl10) < 0) return -1;
	}
	PL = b->gi + i;
	for (i = 0, k = 0; i < 4; ++i)
		for (j = i; j < 4; ++j)
			map[k++] = seq_nt16rev[1<<i|1<<j];
	for (j = 0; j < 3; ++j) // check if ref hom is the most probable in all members
		if (((uint8_t*)PL->data)[j * PL->len] != 0) break;
	if (j < 3) { // we need to go through the complex procedure
		uint8_t *g[3];
		int minc = 1<<30, minc_j = -1, minf = 0, gtf = 0, gtc = 0;
		g[0] = gl10;
		g[1] = gl10 + 10;
		g[2] = gl10 + 20;
		for (j = 1; j <= (int)prep[0]; ++j) { // compute LK with constraint
			int sum = g[0][prep[j]&0xff] + g[1][prep[j]>>8&0xff] + g[2][prep[j]>>16&0xff];
			if (sum < minc) minc = sum, minc_j = j;
		}
		gtc |= map[prep[minc_j]&0xff]; gtc |= map[prep[minc_j]>>8&0xff]<<8; gtc |= map[prep[minc_j]>>16]<<16;
		for (j = 0; j < 3; ++j) { // compute LK without constraint
			int min = 1<<30, min_k = -1;
			for (k = 0; k < 10; ++k)
				if (g[j][k] < min) min = g[j][k], min_k = k;
			gtf |= map[min_k]<<(j*8);
			minf += min;
		}
		*llr = minc - minf; *gt = (int64_t)gtc<<32 | gtf;
	} else *llr = 0, *gt = -1;
	return 0;
}

int bcf_pair_call(const bcf1_t *b)
{
	int i, j, k;
	const bcf_ginfo_t *PL;
	if (b->n_smpl != 2) return -1; // not a pair
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) break;
	if (i == b->n_gi) return -1; // no PL
	PL = b->gi + i;
	for (j = 0; j < 2; ++j) // check if ref hom is the most probable in all members
		if (((uint8_t*)PL->data)[j * PL->len] != 0) break;
	if (j < 2) { // we need to go through the complex procedure
		uint8_t *g[2];
		int minc = 1<<30, minf = 0;
		g[0] = PL->data;
		g[1] = (uint8_t*)PL->data + PL->len;
		for (j = 0; j < PL->len; ++j) // compute LK with constraint
			minc = minc < g[0][j] + g[1][j]? minc : g[0][j] + g[1][j];
		for (j = 0; j < 2; ++j) { // compute LK without constraint
			int min = 1<<30;
			for (k = 0; k < PL->len; ++k)
				min = min < g[j][k]? min : g[j][k];
			minf += min;
		}
		return minc - minf;
	} else return 0;
}

int bcf_min_diff(const bcf1_t *b)
{
	int i, min = 1<<30;
	const bcf_ginfo_t *PL;
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) break;
	if (i == b->n_gi) return -1; // no PL
	PL = b->gi + i;
	for (i = 0; i < b->n_smpl; ++i) {
		int m1, m2, j;
		const uint8_t *p = (uint8_t*)PL->data;
		m1 = m2 = 1<<30;
		for (j = 0; j < PL->len; ++j) {
			if ((int)p[j] < m1) m2 = m1, m1 = p[j];
			else if ((int)p[j] < m2) m2 = p[j];
		}
		min = min < m2 - m1? min : m2 - m1;
	}
	return min;
}
