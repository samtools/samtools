#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bcf.h"

static double g_q2p[256];

#define LD_ITER_MAX 50
#define LD_ITER_EPS 1e-4

#define _G1(h, k) ((h>>1&1) + (k>>1&1))
#define _G2(h, k) ((h&1) + (k&1))

// 0: the previous site; 1: the current site
static int freq_iter(int n, double *pdg[2], double f[4])
{
	double ff[4];
	int i, k, h;
	memset(ff, 0, 4 * sizeof(double));
	for (i = 0; i < n; ++i) {
		double *p[2], sum, tmp;
		p[0] = pdg[0] + i * 3; p[1] = pdg[1] + i * 3;
		for (k = 0, sum = 0.; k < 4; ++k)
			for (h = 0; h < 4; ++h)
				sum += f[k] * f[h] * p[0][_G1(k,h)] * p[1][_G2(k,h)];
		for (k = 0; k < 4; ++k) {
			tmp = f[0] * (p[0][_G1(0,k)] * p[1][_G2(0,k)] + p[0][_G1(k,0)] * p[1][_G2(k,0)])
				+ f[1] * (p[0][_G1(1,k)] * p[1][_G2(1,k)] + p[0][_G1(k,1)] * p[1][_G2(k,1)])
				+ f[2] * (p[0][_G1(2,k)] * p[1][_G2(2,k)] + p[0][_G1(k,2)] * p[1][_G2(k,2)])
				+ f[3] * (p[0][_G1(3,k)] * p[1][_G2(3,k)] + p[0][_G1(k,3)] * p[1][_G2(k,3)]);
			ff[k] += f[k] * tmp / sum;
		}
	}
	for (k = 0; k < 4; ++k) f[k] = ff[k] / (2 * n);
	return 0;
}

double bcf_ld_freq(const bcf1_t *b0, const bcf1_t *b1, double f[4])
{
	const bcf1_t *b[2];
	uint8_t *PL[2];
	int i, j, PL_len[2], n_smpl;
	double *pdg[2], flast[4], r;
	// initialize g_q2p if necessary
	if (g_q2p[0] == 0.)
		for (i = 0; i < 256; ++i)
			g_q2p[i] = pow(10., -i / 10.);
	// initialize others
	if (b0->n_smpl != b1->n_smpl) return -1; // different number of samples
	n_smpl = b0->n_smpl;
	b[0] = b0; b[1] = b1;
	f[0] = f[1] = f[2] = f[3] = -1.;
	if (b[0]->n_alleles < 2 || b[1]->n_alleles < 2) return -1; // one allele only
	// set PL and PL_len
	for (j = 0; j < 2; ++j) {
		const bcf1_t *bj = b[j];
		for (i = 0; i < bj->n_gi; ++i) {
			if (bj->gi[i].fmt == bcf_str2int("PL", 2)) {
				PL[j] = (uint8_t*)bj->gi[i].data;
				PL_len[j] = bj->gi[i].len;
				break;
			}
		}
		if (i == bj->n_gi) return -1; // no PL
	}
	// fill pdg[2]
	pdg[0] = malloc(3 * n_smpl * sizeof(double));
	pdg[1] = malloc(3 * n_smpl * sizeof(double));
	for (j = 0; j < 2; ++j) {
		for (i = 0; i < n_smpl; ++i) {
			const uint8_t *pi = PL[j] + i * PL_len[j];
			double *p = pdg[j] + i * 3;
			p[0] = g_q2p[pi[b[j]->n_alleles]]; p[1] = g_q2p[pi[1]]; p[2] = g_q2p[pi[0]];
		}
	}
	// iteration
	f[0] = f[1] = f[2] = f[3] = 0.25; // this is a really bad guess...
	for (j = 0; j < LD_ITER_MAX; ++j) {
		double eps = 0;
		memcpy(flast, f, 4 * sizeof(double));
		freq_iter(n_smpl, pdg, f);
		for (i = 0; i < 4; ++i) {
			double x = fabs(f[0] - flast[0]);
			if (x > eps) eps = x;
		}
		if (eps < LD_ITER_EPS) break;
	}
	// free
	free(pdg[0]); free(pdg[1]);
	{ // calculate r^2
		double p[2], q[2], D;
		p[0] = f[0] + f[1]; q[0] = 1 - p[0];
		p[1] = f[0] + f[2]; q[1] = 1 - p[1];
		D = f[0] * f[3] - f[1] * f[2];
		r = sqrt(D * D / (p[0] * p[1] * q[0] * q[1]));
		// fprintf(stderr, "R(%lf,%lf,%lf,%lf)=%lf\n", f[0], f[1], f[2], f[3], r2);
		if (isnan(r)) r = -1.;
	}
	return r;
}
