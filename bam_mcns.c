#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "bam_mcns.h"

#define MC_MIN_QUAL 13
#define MC_AVG_ERR 0.007
#define MC_MAX_SUMQ 3000
#define MC_MAX_SUMQP 1e-300
#define MC_MAX_EM_ITER 16
#define MC_EM_EPS 1e-4

struct __mc_aux_t {
	int n, M;
	int ref, alt, alt2;
	double *q2p, *pdg; // pdg -> P(D|g)
	double *phi, *alpha, *CMk; // CMk=\binom{M}{k}
	double *z, *zswap; // aux for afs
	double *afs, *afs1; // afs: accumulative AFS; afs1: site posterior distribution
	int *qsum, *bcnt;
};

static void precal_alpha(mc_aux_t *ma) // \alpha[k]=\binom{M}{k}\sum_l\phi_l(l/M)^k(1-l/M)^{M-k}
{
	int k, l;
	memset(ma->alpha, 0, sizeof(double) * (ma->M + 1));
	for (l = 0; l <= ma->M; ++l)
		ma->alpha[0] += ma->phi[l] * pow(1. - (double)l / ma->M, ma->M);
	for (k = 1; k < ma->M; ++k) {
		for (l = 1; l < ma->M; ++l) { // for k!=0 and k!=ma->M, l=0 and l=ma->M leads to zero
			double x = exp((lgamma(ma->M + 1) - lgamma(k + 1) - lgamma(ma->M - k + 1))
						   + k * log((double)l / ma->M)
						   + (ma->M - k) * log(1. - (double)l / ma->M));
			ma->alpha[k] += x * ma->phi[l];
		}
	}
	for (l = 0; l <= ma->M; ++l)
		ma->alpha[ma->M] += ma->phi[l] * pow((double)l / ma->M, ma->M);
	fflush(stdout);
}

void mc_init_prior(mc_aux_t *ma, int type, double theta)
{
	int i;
	if (type == MC_PTYPE_COND2) {
		for (i = 0; i <= 2 * ma->n; ++i)
			ma->phi[i] = 2. * (i + 1) / (2 * ma->n + 1) / (2 * ma->n + 2);
	} else if (type == MC_PTYPE_FLAT) {
		for (i = 0; i <= ma->M; ++i)
			ma->phi[i] = 1. / (ma->M + 1);
	} else {
		double sum;
		for (i = 0, sum = 0.; i < 2 * ma->n; ++i)
			sum += (ma->phi[i] = theta / (2 * ma->n - i));
		ma->phi[2 * ma->n] = 1. - sum;
	}
	precal_alpha(ma);
}

mc_aux_t *mc_init(int n) // FIXME: assuming diploid
{
	mc_aux_t *ma;
	int i;
	ma = calloc(1, sizeof(mc_aux_t));
	ma->n = n; ma->M = 2 * n;
	ma->q2p = calloc(MC_MAX_SUMQ + 1, sizeof(double));
	ma->qsum = calloc(4 * ma->n, sizeof(int));
	ma->bcnt = calloc(4 * ma->n, sizeof(int));
	ma->pdg = calloc(3 * ma->n, sizeof(double));
	ma->phi = calloc(ma->M + 1, sizeof(double));
	ma->alpha = calloc(ma->M + 1, sizeof(double));
	ma->CMk = calloc(ma->M + 1, sizeof(double));
	ma->z = calloc(2 * ma->n + 1, sizeof(double));
	ma->zswap = calloc(2 * ma->n + 1, sizeof(double));
	ma->afs = calloc(2 * ma->n + 1, sizeof(double));
	ma->afs1 = calloc(2 * ma->n + 1, sizeof(double));
	for (i = 0; i <= MC_MAX_SUMQ; ++i)
		ma->q2p[i] = pow(10., -i / 10.);
	for (i = 0; i <= ma->M; ++i)
		ma->CMk[i] = exp(lgamma(ma->M + 1) - lgamma(i + 1) - lgamma(ma->M - i + 1));
	mc_init_prior(ma, MC_PTYPE_FULL, 1e-3); // the simplest prior
	return ma;
}

void mc_destroy(mc_aux_t *ma)
{
	if (ma) {
		free(ma->qsum); free(ma->bcnt);
		free(ma->q2p); free(ma->pdg);
		free(ma->phi); free(ma->alpha); free(ma->CMk);
		free(ma->z); free(ma->zswap);
		free(ma->afs); free(ma->afs1);
		free(ma);
	}
}

static int sum_err(int *n, const bam_pileup1_t **plp, mc_aux_t *ma)
{
	int i, j, tot = 0;
	memset(ma->qsum, 0, sizeof(int) * 4 * ma->n);
	memset(ma->bcnt, 0, sizeof(int) * 4 * ma->n);
	for (j = 0; j < ma->n; ++j) {
		int *qsum = ma->qsum + j * 4;
		int *bcnt = ma->bcnt + j * 4;
		for (i = 0; i < n[j]; ++i) {
			const bam_pileup1_t *p = plp[j] + i;
			int q, b;
			if (p->is_del || (p->b->core.flag&BAM_FUNMAP)) continue;
			q = bam1_qual(p->b)[p->qpos];
			if (p->b->core.qual < q) q = p->b->core.qual;
			if (q < MC_MIN_QUAL) continue; // small qual
			b = bam_nt16_nt4_table[(int)bam1_seqi(bam1_seq(p->b), p->qpos)];
			if (b > 3) continue; // N
			qsum[b] += q;
			++bcnt[b];
			++tot;
		}
	}
	return tot;
}

static void set_allele(int ref, mc_aux_t *ma)
{
	int i, j, sum[4], tmp;
	sum[0] = sum[1] = sum[2] = sum[3] = 0;
	for (i = 0; i < ma->n; ++i)
		for (j = 0; j < 4; ++j)
			sum[j] += ma->qsum[i * 4 + j];
	for (j = 0; j < 4; ++j) sum[j] = sum[j]<<2 | j;
	for (i = 1; i < 4; ++i) // insertion sort
		for (j = i; j > 0 && sum[j] < sum[j-1]; --j)
			tmp = sum[j], sum[j] = sum[j-1], sum[j-1] = tmp;
	ma->ref = sum[3]&3; ma->alt = sum[2]&3; ma->alt2 = -1;
	if (ma->ref != ref) { // the best base is not ref
		if (ref >= 0 && ref <= 3) { // ref is not N
			if (ma->alt == ref) tmp = ma->ref, ma->ref = ma->alt, ma->alt = tmp; // then switch alt and ref
			else ma->alt2 = ma->alt, ma->alt = ma->ref, ma->ref = ref; // then set ref as ref
		} else ma->alt2 = ma->alt, ma->alt = ma->ref, ma->ref = sum[0]&3; // then set the weakest as ref
	}
}

static void cal_pdg(mc_aux_t *ma)
{
	int i, j;
	for (j = 0; j < ma->n; ++j) {
		int pi[3], *qsum, *bcnt;
		double *pdg = ma->pdg + j * 3;
		qsum = ma->qsum + j * 4;
		bcnt = ma->bcnt + j * 4;
		pi[1] = 3 * (bcnt[ma->ref] + bcnt[ma->alt]);
		pi[0] = qsum[ma->ref];
		pi[2] = qsum[ma->alt];
		for (i = 0; i < 3; ++i)
			pdg[i] = pi[i] > MC_MAX_SUMQ? MC_MAX_SUMQP : ma->q2p[pi[i]];
	}
}
// this calculates the naive allele frequency and Nielsen's frequency
static double mc_freq0(const mc_aux_t *ma, double *_f)
{
	int i, cnt;
	double f, f_nielsen, w_sum;
	*_f = -1.;
	for (i = cnt = 0, f = f_nielsen = w_sum = 0.; i < ma->n; ++i) {
		int *bcnt = ma->bcnt + i * 4;
		int x = bcnt[ma->ref] + bcnt[ma->alt];
		if (x) {
			double w, p;
			++cnt;
			f += (double)bcnt[ma->ref] / x;
			p = (bcnt[ma->ref] - MC_AVG_ERR * x) / (1. - 2. * MC_AVG_ERR) / x;
			w = 2. * x / (1. + x);
			w_sum += w;
			f_nielsen += p * w;
		}
	}
	if (cnt) {
		f_nielsen /= w_sum;
		if (f_nielsen < 0.) f_nielsen = 0.;
		if (f_nielsen > 1.) f_nielsen = 1.;
		*_f = f_nielsen;
		return f / cnt;
	} else return -1.;
}
// f0 is the reference allele frequency
static double mc_freq_iter(double f0, const mc_aux_t *ma)
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

int mc_call_gt(const mc_aux_t *ma, double f0, int k)
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

static void mc_cal_z(mc_aux_t *ma)
{
	double *z[2], *tmp, *pdg;
	int i, j;
	z[0] = ma->z;
	z[1] = ma->zswap;
	pdg = ma->pdg;
	z[0][0] = 1.; z[0][1] = z[0][2] = 0.;
	for (j = 0; j < ma->n; ++j) {
		int max = (j + 1) * 2;
		double p[3];
		pdg = ma->pdg + j * 3;
		p[0] = pdg[0]; p[1] = 2. * pdg[1]; p[2] = pdg[2];
		z[1][0] = p[0] * z[0][0];
		z[1][1] = p[0] * z[0][1] + p[1] * z[0][0];
		for (i = 2; i <= max; ++i)
			z[1][i] = p[0] * z[0][i] + p[1] * z[0][i-1] + p[2] * z[0][i-2];
		if (j < ma->n - 1) z[1][max+1] = z[1][max+2] = 0.;
//		int k; for (k = 0; k <= max; ++k) printf("%d:%.3lg ", k, z[1][k]); putchar('\n');
		tmp = z[0]; z[0] = z[1]; z[1] = tmp;
	}
	if (z[0] != ma->z) memcpy(ma->z, z[0], sizeof(double) * (2 * ma->n + 1));
}

static double mc_add_afs(mc_aux_t *ma)
{
	int k;
	long double sum = 0.;
	memset(ma->afs1, 0, sizeof(double) * (ma->M + 1));
	mc_cal_z(ma);
	for (k = 0, sum = 0.; k <= ma->M; ++k)
		sum += (long double)ma->alpha[k] * ma->z[k] / ma->CMk[k];
	for (k = 0; k <= ma->M; ++k) {
		ma->afs1[k] = ma->alpha[k] * ma->z[k] / ma->CMk[k] / sum;
		if (isnan(ma->afs1[k]) || isinf(ma->afs1[k])) return -1.;
	}
	for (k = 0, sum = 0.; k <= ma->M; ++k) {
		ma->afs[k] += ma->afs1[k];
		sum += k * ma->afs1[k];
	}
	return sum / ma->M;
}

int mc_cal(int ref, int *n, const bam_pileup1_t **plp, mc_aux_t *ma, mc_rst_t *rst, int level)
{
	int i, tot;
	memset(rst, 0, sizeof(mc_rst_t));
	rst->f_em = rst->f_exp = -1.; rst->ref = rst->alt = -1;
	// precalculation
	tot = sum_err(n, plp, ma);
	if (tot == 0) return 0; // no good bases
	set_allele(ref, ma);
	cal_pdg(ma);
	// set ref/major allele
	rst->ref = ma->ref; rst->alt = ma->alt; rst->alt2 = ma->alt2;
	// calculate naive and Nielsen's freq
	rst->f_naive = mc_freq0(ma, &rst->f_nielsen);
	{ // calculate f_em
		double flast = rst->f_naive;
		for (i = 0; i < MC_MAX_EM_ITER; ++i) {
			rst->f_em = mc_freq_iter(flast, ma);
			if (fabs(rst->f_em - flast) < MC_EM_EPS) break;
			flast = rst->f_em;
		}
	}
	if (level >= 2) {
		rst->f_exp = mc_add_afs(ma);
		rst->p_ref = ma->afs1[ma->M];
	}
	return tot;
}

void mc_dump_afs(mc_aux_t *ma)
{
	int k;
	fprintf(stderr, "[afs]");
	for (k = 0; k <= ma->M; ++k)
		fprintf(stderr, " %d:%.3lf", k, ma->afs[ma->M - k]);
	fprintf(stderr, "\n");
	memset(ma->afs, 0, sizeof(double) * (ma->M + 1));
}
