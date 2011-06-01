#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bcf.h"
#include "kmin.h"

static double g_q2p[256];

#define ITER_MAX 50
#define ITER_TRY 10
#define EPS 1e-5

extern double kf_gammaq(double, double);

/*
	Generic routines
 */
// get the 3 genotype likelihoods
static double *get_pdg3(const bcf1_t *b)
{
	double *pdg;
	const uint8_t *PL = 0;
	int i, PL_len = 0;
	// initialize g_q2p if necessary
	if (g_q2p[0] == 0.)
		for (i = 0; i < 256; ++i)
			g_q2p[i] = pow(10., -i / 10.);
	// set PL and PL_len
	for (i = 0; i < b->n_gi; ++i) {
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
			PL = (const uint8_t*)b->gi[i].data;
			PL_len = b->gi[i].len;
			break;
		}
	}
	if (i == b->n_gi) return 0; // no PL
	// fill pdg
	pdg = malloc(3 * b->n_smpl * sizeof(double));
	for (i = 0; i < b->n_smpl; ++i) {
		const uint8_t *pi = PL + i * PL_len;
		double *p = pdg + i * 3;
		p[0] = g_q2p[pi[2]]; p[1] = g_q2p[pi[1]]; p[2] = g_q2p[pi[0]];
	}
	return pdg;
}

// estimate site allele frequency in a very naive and inaccurate way
static double est_freq(int n, const double *pdg)
{
	int i, gcnt[3], tmp1;
	// get a rough estimate of the genotype frequency
	gcnt[0] = gcnt[1] = gcnt[2] = 0;
	for (i = 0; i < n; ++i) {
		const double *p = pdg + i * 3;
		if (p[0] != 1. || p[1] != 1. || p[2] != 1.) {
			int which = p[0] > p[1]? 0 : 1;
			which = p[which] > p[2]? which : 2;
			++gcnt[which];
		}
	}
	tmp1 = gcnt[0] + gcnt[1] + gcnt[2];
	return (tmp1 == 0)? -1.0 : (.5 * gcnt[1] + gcnt[2]) / tmp1;
}

/*
	Single-locus EM
 */

typedef struct {
	int beg, end;
	const double *pdg;
} minaux1_t;

static double prob1(double f, void *data)
{
	minaux1_t *a = (minaux1_t*)data;
	double p = 1., l = 0., f3[3];
	int i;
//	printf("brent %lg\n", f);
	if (f < 0 || f > 1) return 1e300;
	f3[0] = (1.-f)*(1.-f); f3[1] = 2.*f*(1.-f); f3[2] = f*f;
	for (i = a->beg; i < a->end; ++i) {
		const double *pdg = a->pdg + i * 3;
		p *= pdg[0] * f3[0] + pdg[1] * f3[1] + pdg[2] * f3[2];
		if (p < 1e-200) l -= log(p), p = 1.;
	}
	return l - log(p);
}

// one EM iteration for allele frequency estimate
static double freq_iter(double *f, const double *_pdg, int beg, int end)
{
	double f0 = *f, f3[3], err;
	int i;
//	printf("em %lg\n", *f);
	f3[0] = (1.-f0)*(1.-f0); f3[1] = 2.*f0*(1.-f0); f3[2] = f0*f0;
	for (i = beg, f0 = 0.; i < end; ++i) {
		const double *pdg = _pdg + i * 3;
		f0 += (pdg[1] * f3[1] + 2. * pdg[2] * f3[2])
			/ (pdg[0] * f3[0] + pdg[1] * f3[1] + pdg[2] * f3[2]);
	}
	f0 /= (end - beg) * 2;
	err = fabs(f0 - *f);
	*f = f0;
	return err;
}

/* The following function combines EM and Brent's method. When the signal from
 * the data is strong, EM is faster but sometimes, EM may converge very slowly.
 * When this happens, we switch to Brent's method. The idea is learned from
 * Rasmus Nielsen.
 */
static double freqml(double f0, int beg, int end, const double *pdg)
{
	int i;
	double f;
	for (i = 0, f = f0; i < ITER_TRY; ++i)
		if (freq_iter(&f, pdg, beg, end) < EPS) break;
	if (i == ITER_TRY) { // haven't converged yet; try Brent's method
		minaux1_t a;
		a.beg = beg; a.end = end; a.pdg = pdg;
		kmin_brent(prob1, f0 == f? .5*f0 : f0, f, (void*)&a, EPS, &f);
	}
	return f;
}

// one EM iteration for genotype frequency estimate
static double g3_iter(double g[3], const double *_pdg, int beg, int end)
{
	double err, gg[3];
	int i;
	gg[0] = gg[1] = gg[2] = 0.;
//	printf("%lg,%lg,%lg\n", g[0], g[1], g[2]);
	for (i = beg; i < end; ++i) {
		double sum, tmp[3];
		const double *pdg = _pdg + i * 3;
		tmp[0] = pdg[0] * g[0]; tmp[1] = pdg[1] * g[1]; tmp[2] = pdg[2] * g[2];
		sum = (tmp[0] + tmp[1] + tmp[2]) * (end - beg);
		gg[0] += tmp[0] / sum; gg[1] += tmp[1] / sum; gg[2] += tmp[2] / sum;
	}
	err = fabs(gg[0] - g[0]) > fabs(gg[1] - g[1])? fabs(gg[0] - g[0]) : fabs(gg[1] - g[1]);
	err = err > fabs(gg[2] - g[2])? err : fabs(gg[2] - g[2]);
	g[0] = gg[0]; g[1] = gg[1]; g[2] = gg[2];
	return err;
}

// perform likelihood ratio test
static double lk_ratio_test(int n, int n1, const double *pdg, double f3[3][3])
{
	double r;
	int i;
	for (i = 0, r = 1.; i < n1; ++i) {
		const double *p = pdg + i * 3;
		r *= (p[0] * f3[1][0] + p[1] * f3[1][1] + p[2] * f3[1][2])
			/ (p[0] * f3[0][0] + p[1] * f3[0][1] + p[2] * f3[0][2]);
	}
	for (; i < n; ++i) {
		const double *p = pdg + i * 3;
		r *= (p[0] * f3[2][0] + p[1] * f3[2][1] + p[2] * f3[2][2])
			/ (p[0] * f3[0][0] + p[1] * f3[0][1] + p[2] * f3[0][2]);
	}
	return r;
}

// x[0]: ref frequency
// x[1..3]: alt-alt, alt-ref, ref-ref frequenc
// x[4]: HWE P-value
// x[5..6]: group1 freq, group2 freq
// x[7]: 1-degree P-value
// x[8]: 2-degree P-value
// x[9]: 1-degree P-value with freq estimated from genotypes
int bcf_em1(const bcf1_t *b, int n1, int flag, double x[10])
{
	double *pdg;
	int i, n, n2;
	if (b->n_alleles < 2) return -1; // one allele only
	// initialization
	if (n1 < 0 || n1 > b->n_smpl) n1 = 0;
	if (flag & 1<<7) flag |= 7<<5; // compute group freq if LRT is required
	if (flag & 0xf<<1) flag |= 0xf<<1;
	n = b->n_smpl; n2 = n - n1;
	pdg = get_pdg3(b);
	if (pdg == 0) return -1;
	for (i = 0; i < 10; ++i) x[i] = -1.; // set to negative
	{
		if ((x[0] = est_freq(n, pdg)) < 0.) {
			free(pdg);
			return -1; // no data
		}
		x[0] = freqml(x[0], 0, n, pdg);
	}
	if (flag & (0xf<<1|3<<8)) { // estimate the genotype frequency and test HWE
		double *g = x + 1, f3[3], r;
		f3[0] = g[0] = (1 - x[0]) * (1 - x[0]);
		f3[1] = g[1] = 2 * x[0] * (1 - x[0]);
		f3[2] = g[2] = x[0] * x[0];
		for (i = 0; i < ITER_MAX; ++i)
			if (g3_iter(g, pdg, 0, n) < EPS) break;
		// Hardy-Weinberg equilibrium (HWE)
		for (i = 0, r = 1.; i < n; ++i) {
			double *p = pdg + i * 3;
			r *= (p[0] * g[0] + p[1] * g[1] + p[2] * g[2]) / (p[0] * f3[0] + p[1] * f3[1] + p[2] * f3[2]);
		}
		x[4] = kf_gammaq(.5, log(r));
	}
	if ((flag & 7<<5) && n1 > 0 && n1 < n) { // group frequency
		x[5] = freqml(x[0], 0, n1, pdg);
		x[6] = freqml(x[0], n1, n, pdg);
	}
	if ((flag & 1<<7) && n1 > 0 && n1 < n) { // 1-degree P-value
		double f[3], f3[3][3];
		f[0] = x[0]; f[1] = x[5]; f[2] = x[6];
		for (i = 0; i < 3; ++i)
			f3[i][0] = (1-f[i])*(1-f[i]), f3[i][1] = 2*f[i]*(1-f[i]), f3[i][2] = f[i]*f[i];
		x[7] = kf_gammaq(.5, log(lk_ratio_test(n, n1, pdg, f3)));
	}
	if ((flag & 3<<8) && n1 > 0 && n1 < n) { // 2-degree P-value
		double g[3][3], tmp;
		for (i = 0; i < 3; ++i) memcpy(g[i], x + 1, 3 * sizeof(double));
		for (i = 0; i < ITER_MAX; ++i)
			if (g3_iter(g[1], pdg, 0, n1) < EPS) break;
		for (i = 0; i < ITER_MAX; ++i)
			if (g3_iter(g[2], pdg, n1, n) < EPS) break;
		tmp = log(lk_ratio_test(n, n1, pdg, g));
		x[8] = kf_gammaq(1., tmp);
		x[9] = kf_gammaq(.5, tmp);
	}
	// free
	free(pdg);
	return 0;
}

/*
	Two-locus EM (LD)
 */

#define _G1(h, k) ((h>>1&1) + (k>>1&1))
#define _G2(h, k) ((h&1) + (k&1))

// 0: the previous site; 1: the current site
static int pair_freq_iter(int n, double *pdg[2], double f[4])
{
	double ff[4];
	int i, k, h;
//	printf("%lf,%lf,%lf,%lf\n", f[0], f[1], f[2], f[3]);
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

double bcf_pair_freq(const bcf1_t *b0, const bcf1_t *b1, double f[4])
{
	const bcf1_t *b[2];
	int i, j, n_smpl;
	double *pdg[2], flast[4], r, f0[2];
	// initialize others
	if (b0->n_smpl != b1->n_smpl) return -1; // different number of samples
	n_smpl = b0->n_smpl;
	b[0] = b0; b[1] = b1;
	f[0] = f[1] = f[2] = f[3] = -1.;
	if (b[0]->n_alleles < 2 || b[1]->n_alleles < 2) return -1; // one allele only
	pdg[0] = get_pdg3(b0); pdg[1] = get_pdg3(b1);
	if (pdg[0] == 0 || pdg[1] == 0) {
		free(pdg[0]); free(pdg[1]);
		return -1;
	}
	// set the initial value
	f0[0] = est_freq(n_smpl, pdg[0]);
	f0[1] = est_freq(n_smpl, pdg[1]);
	f[0] = (1 - f0[0]) * (1 - f0[1]); f[3] = f0[0] * f0[1];
	f[1] = (1 - f0[0]) * f0[1]; f[2] = f0[0] * (1 - f0[1]);
	// iteration
	for (j = 0; j < ITER_MAX; ++j) {
		double eps = 0;
		memcpy(flast, f, 4 * sizeof(double));
		pair_freq_iter(n_smpl, pdg, f);
		for (i = 0; i < 4; ++i) {
			double x = fabs(f[i] - flast[i]);
			if (x > eps) eps = x;
		}
		if (eps < EPS) break;
	}
	// free
	free(pdg[0]); free(pdg[1]);
	{ // calculate r^2
		double p[2], q[2], D;
		p[0] = f[0] + f[1]; q[0] = 1 - p[0];
		p[1] = f[0] + f[2]; q[1] = 1 - p[1];
		D = f[0] * f[3] - f[1] * f[2];
		r = sqrt(D * D / (p[0] * p[1] * q[0] * q[1]));
//		printf("R(%lf,%lf,%lf,%lf)=%lf\n", f[0], f[1], f[2], f[3], r);
		if (isnan(r)) r = -1.;
	}
	return r;
}
