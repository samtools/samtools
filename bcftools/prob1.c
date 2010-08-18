#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "prob1.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define MC_AVG_ERR 0.007
#define MC_MAX_EM_ITER 16
#define MC_EM_EPS 1e-4

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

struct __bcf_p1aux_t {
	int n, M;
	double *q2p, *pdg; // pdg -> P(D|g)
	double *phi;
	double *z, *zswap; // aux for afs
	double *afs, *afs1; // afs: accumulative AFS; afs1: site posterior distribution
	const uint8_t *PL; // point to PL
	int PL_len;
};

void bcf_p1_init_prior(bcf_p1aux_t *ma, int type, double theta)
{
	int i;
	if (type == MC_PTYPE_COND2) {
		for (i = 0; i <= ma->M; ++i)
			ma->phi[i] = 2. * (i + 1) / (ma->M + 1) / (ma->M + 2);
	} else if (type == MC_PTYPE_FLAT) {
		for (i = 0; i <= ma->M; ++i)
			ma->phi[i] = 1. / (ma->M + 1);
	} else {
		double sum;
		for (i = 0, sum = 0.; i < ma->M; ++i)
			sum += (ma->phi[i] = theta / (ma->M - i));
		ma->phi[ma->M] = 1. - sum;
	}
}

int bcf_p1_read_prior(bcf_p1aux_t *ma, const char *fn)
{
	gzFile fp;
	kstring_t s;
	kstream_t *ks;
	long double sum;
	int dret, k;
	memset(&s, 0, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	memset(ma->phi, 0, sizeof(double) * (ma->M + 1));
	while (ks_getuntil(ks, '\n', &s, &dret) >= 0) {
		if (strstr(s.s, "[afs] ") == s.s) {
			char *p = s.s + 6;
			for (k = 0; k <= ma->M; ++k) {
				int x;
				double y;
				x = strtol(p, &p, 10);
				if (x != k && (errno == EINVAL || errno == ERANGE)) return -1;
				++p;
				y = strtod(p, &p);
				if (y == 0. && (errno == EINVAL || errno == ERANGE)) return -1;
				ma->phi[ma->M - k] += y;
			}
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	free(s.s);
	for (sum = 0., k = 0; k <= ma->M; ++k) sum += ma->phi[k];
	fprintf(stderr, "[prior]");
	for (k = 0; k <= ma->M; ++k) ma->phi[k] /= sum;
	for (k = 0; k <= ma->M; ++k) fprintf(stderr, " %d:%.3lg", k, ma->phi[ma->M - k]);
	fputc('\n', stderr);
	return 0;
}

bcf_p1aux_t *bcf_p1_init(int n) // FIXME: assuming diploid
{
	bcf_p1aux_t *ma;
	int i;
	ma = calloc(1, sizeof(bcf_p1aux_t));
	ma->n = n; ma->M = 2 * n;
	ma->q2p = calloc(256, sizeof(double));
	ma->pdg = calloc(3 * ma->n, sizeof(double));
	ma->phi = calloc(ma->M + 1, sizeof(double));
	ma->z = calloc(2 * ma->n + 1, sizeof(double));
	ma->zswap = calloc(2 * ma->n + 1, sizeof(double));
	ma->afs = calloc(2 * ma->n + 1, sizeof(double));
	ma->afs1 = calloc(2 * ma->n + 1, sizeof(double));
	for (i = 0; i < 256; ++i)
		ma->q2p[i] = pow(10., -i / 10.);
	bcf_p1_init_prior(ma, MC_PTYPE_FULL, 1e-3); // the simplest prior
	return ma;
}

void bcf_p1_destroy(bcf_p1aux_t *ma)
{
	if (ma) {
		free(ma->q2p); free(ma->pdg);
		free(ma->phi);
		free(ma->z); free(ma->zswap);
		free(ma->afs); free(ma->afs1);
		free(ma);
	}
}

#define char2int(s) (((int)s[0])<<8|s[1])

static int cal_pdg(const bcf1_t *b, bcf_p1aux_t *ma)
{
	int i, j, k;
	long *p, tmp;
	p = alloca(b->n_alleles * sizeof(long));
	memset(p, 0, sizeof(long) * b->n_alleles);
	for (j = 0; j < ma->n; ++j) {
		const uint8_t *pi = ma->PL + j * ma->PL_len;
		double *pdg = ma->pdg + j * 3;
		pdg[0] = ma->q2p[pi[b->n_alleles]]; pdg[1] = ma->q2p[pi[1]]; pdg[2] = ma->q2p[pi[0]];
		for (i = k = 0; i < b->n_alleles; ++i) {
			p[i] += (int)pi[k];
			k += b->n_alleles - i;
		}
	}
	for (i = 0; i < b->n_alleles; ++i) p[i] = p[i]<<4 | i;
	for (i = 1; i < b->n_alleles; ++i) // insertion sort
		for (j = i; j > 0 && p[j] < p[j-1]; --j)
			tmp = p[j], p[j] = p[j-1], p[j-1] = tmp;
	for (i = b->n_alleles - 1; i >= 0; --i)
		if ((p[i]&0xf) == 0) break;
	return i;
}
// f0 is the reference allele frequency
static double mc_freq_iter(double f0, const bcf_p1aux_t *ma)
{
	double f, f3[3];
	int i;
	f3[0] = (1.-f0)*(1.-f0); f3[1] = 2.*f0*(1.-f0); f3[2] = f0*f0;
	for (i = 0, f = 0.; i < ma->n; ++i) {
		double *pdg;
		pdg = ma->pdg + i * 3;
		f += (pdg[1] * f3[1] + 2. * pdg[2] * f3[2])
			/ (pdg[0] * f3[0] + pdg[1] * f3[1] + pdg[2] * f3[2]);
	}
	f /= ma->n * 2.;
	return f;
}

int bcf_p1_call_gt(const bcf_p1aux_t *ma, double f0, int k)
{
	double sum, g[3];
	double max, f3[3], *pdg = ma->pdg + k * 3;
	int q, i, max_i;
	f3[0] = (1.-f0)*(1.-f0); f3[1] = 2.*f0*(1.-f0); f3[2] = f0*f0;
	for (i = 0, sum = 0.; i < 3; ++i)
		sum += (g[i] = pdg[i] * f3[i]);
	for (i = 0, max = -1., max_i = 0; i < 3; ++i) {
		g[i] /= sum;
		if (g[i] > max) max = g[i], max_i = i;
	}
	max = 1. - max;
	if (max < 1e-308) max = 1e-308;
	q = (int)(-3.434 * log(max) + .499);
	if (q > 99) q = 99;
	return q<<2|max_i;
}

#define TINY 1e-20

static void mc_cal_y(bcf_p1aux_t *ma)
{
	double *z[2], *tmp, *pdg;
	int k, j, last_min, last_max;
	z[0] = ma->z;
	z[1] = ma->zswap;
	pdg = ma->pdg;
	memset(z[0], 0, sizeof(double) * (ma->M + 1));
	memset(z[1], 0, sizeof(double) * (ma->M + 1));
	z[0][0] = 1.;
	last_min = last_max = 0;
	for (j = 0; j < ma->n; ++j) {
		int _min = last_min, _max = last_max;
		double p[3], sum;
		pdg = ma->pdg + j * 3;
		p[0] = pdg[0]; p[1] = 2. * pdg[1]; p[2] = pdg[2];
		for (; _min < _max && z[0][_min] < TINY; ++_min) z[0][_min] = z[1][_min] = 0.;
		for (; _max > _min && z[0][_max] < TINY; --_max) z[0][_max] = z[1][_max] = 0.;
		_max += 2;
		if (_min == 0) 
			k = 0, z[1][k] = (2*j+2-k)*(2*j-k+1) * p[0] * z[0][k];
		if (_min <= 1)
			k = 1, z[1][k] = (2*j+2-k)*(2*j-k+1) * p[0] * z[0][k] + k*(2*j+2-k) * p[1] * z[0][k-1];
		for (k = _min < 2? 2 : _min; k <= _max; ++k)
			z[1][k] = (2*j+2-k)*(2*j-k+1) * p[0] * z[0][k]
				+ k*(2*j+2-k) * p[1] * z[0][k-1]
				+ k*(k-1)* p[2] * z[0][k-2];
		for (k = _min, sum = 0.; k <= _max; ++k) sum += z[1][k];
		for (k = _min; k <= _max; ++k) z[1][k] /= sum;
		if (_min >= 1) z[1][_min-1] = 0.;
		if (_min >= 2) z[1][_min-2] = 0.;
		if (j < ma->n - 1) z[1][_max+1] = z[1][_max+2] = 0.;
		tmp = z[0]; z[0] = z[1]; z[1] = tmp;
		last_min = _min; last_max = _max;
	}
	if (z[0] != ma->z) memcpy(ma->z, z[0], sizeof(double) * (ma->M + 1));
}

static double mc_cal_afs(bcf_p1aux_t *ma)
{
	int k;
	long double sum = 0.;
	memset(ma->afs1, 0, sizeof(double) * (ma->M + 1));
	mc_cal_y(ma);
	for (k = 0, sum = 0.; k <= ma->M; ++k)
		sum += (long double)ma->phi[k] * ma->z[k];
	for (k = 0; k <= ma->M; ++k) {
		ma->afs1[k] = ma->phi[k] * ma->z[k] / sum;
		if (isnan(ma->afs1[k]) || isinf(ma->afs1[k])) return -1.;
	}
	for (k = 0, sum = 0.; k <= ma->M; ++k) {
		ma->afs[k] += ma->afs1[k];
		sum += k * ma->afs1[k];
	}
	return sum / ma->M;
}

static long double p1_cal_g3(bcf_p1aux_t *p1a, double g[3])
{
	long double pd = 0., g2[3];
	int i, k;
	memset(g2, 0, sizeof(long double) * 3);
	for (k = 0; k < p1a->M; ++k) {
		double f = (double)k / p1a->M, f3[3], g1[3];
		long double z = 1.;
		g1[0] = g1[1] = g1[2] = 0.;
		f3[0] = (1. - f) * (1. - f); f3[1] = 2. * f * (1. - f); f3[2] = f * f;
		for (i = 0; i < p1a->n; ++i) {
			double *pdg = p1a->pdg + i * 3;
			double x = pdg[0] * f3[0] + pdg[1] * f3[1] + pdg[2] * f3[2];
			z *= x;
			g1[0] += pdg[0] * f3[0] / x;
			g1[1] += pdg[1] * f3[1] / x;
			g1[2] += pdg[2] * f3[2] / x;
		}
		pd += p1a->phi[k] * z;
		for (i = 0; i < 3; ++i)
			g2[i] += p1a->phi[k] * z * g1[i];
	}
	for (i = 0; i < 3; ++i) g[i] = g2[i] / pd;
	return pd;
}

int bcf_p1_cal(bcf1_t *b, bcf_p1aux_t *ma, bcf_p1rst_t *rst)
{
	int i, k;
	long double sum = 0.;
	// set PL and PL_len
	for (i = 0; i < b->n_gi; ++i) {
		if (b->gi[i].fmt == char2int("PL")) {
			ma->PL = (uint8_t*)b->gi[i].data;
			ma->PL_len = b->gi[i].len;
			break;
		}
	}
	if (b->n_alleles < 2) return -1; // FIXME: find a better solution
	// 
	rst->rank0 = cal_pdg(b, ma);
	rst->f_exp = mc_cal_afs(ma);
	rst->p_ref = ma->afs1[ma->M];
	// calculate f_flat and f_em
	for (k = 0, sum = 0.; k <= ma->M; ++k)
		sum += (long double)ma->z[k];
	rst->f_flat = 0.;
	for (k = 0; k <= ma->M; ++k) {
		double p = ma->z[k] / sum;
		rst->f_flat += k * p;
	}
	rst->f_flat /= ma->M;
	{ // calculate f_em
		double flast = rst->f_flat;
		for (i = 0; i < MC_MAX_EM_ITER; ++i) {
			rst->f_em = mc_freq_iter(flast, ma);
			if (fabs(rst->f_em - flast) < MC_EM_EPS) break;
			flast = rst->f_em;
		}
	}
	p1_cal_g3(ma, rst->g);
	return 0;
}

void bcf_p1_dump_afs(bcf_p1aux_t *ma)
{
	int k;
	fprintf(stderr, "[afs]");
	for (k = 0; k <= ma->M; ++k)
		fprintf(stderr, " %d:%.3lf", k, ma->afs[ma->M - k]);
	fprintf(stderr, "\n");
	memset(ma->afs, 0, sizeof(double) * (ma->M + 1));
}
