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
	double *alpha, *beta;
	double *z, *zswap; // aux for afs
	double *afs, *afs1; // afs: accumulative AFS; afs1: site posterior distribution
	int *qsum, *bcnt;
};

void mc_init_prior(mc_aux_t *ma, int type, double theta)
{
	int i;
	if (type == MC_PTYPE_COND2) {
		for (i = 0; i <= 2 * ma->n; ++i)
			ma->alpha[i] = 2. * (i + 1) / (2 * ma->n + 1) / (2 * ma->n + 2);
	} else if (type == MC_PTYPE_FLAT) {
		for (i = 0; i <= ma->M; ++i)
			ma->alpha[i] = 1. / (ma->M + 1);
	} else {
		double sum;
		for (i = 0, sum = 0.; i < 2 * ma->n; ++i)
			sum += (ma->alpha[i] = theta / (2 * ma->n - i));
		ma->alpha[2 * ma->n] = 1. - sum;
	}
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
	ma->alpha = calloc(2 * ma->n + 1, sizeof(double));
	ma->beta = calloc((2 * ma->n + 1) * 3, sizeof(double));
	ma->z = calloc(2 * ma->n + 1, sizeof(double));
	ma->zswap = calloc(2 * ma->n + 1, sizeof(double));
	ma->afs = calloc(2 * ma->n + 1, sizeof(double));
	ma->afs1 = calloc(2 * ma->n + 1, sizeof(double));
	for (i = 0; i <= MC_MAX_SUMQ; ++i)
		ma->q2p[i] = pow(10., -i / 10.);
	for (i = 0; i <= ma->M; ++i) { // beta[k][g]=P(g|k/M)
		double *bi = ma->beta + 3 * i;
		double f = (double)i / ma->M;
		bi[0] = (1. - f) * (1. - f);
		bi[1] = 2 * f * (1. - f);
		bi[2] = f * f;
	}
	mc_init_prior(ma, MC_PTYPE_FULL, 1e-3); // the simplest prior
	return ma;
}

void mc_destroy(mc_aux_t *ma)
{
	if (ma) {
		free(ma->qsum); free(ma->bcnt);
		free(ma->q2p); free(ma->pdg);
		free(ma->alpha); free(ma->beta);
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
	if (ma->ref != ref) {
		if (ma->alt == ref) tmp = ma->ref, ma->ref = ma->alt, ma->alt = tmp;
		else ma->alt2 = ma->alt, ma->alt = ma->ref, ma->ref = ref;
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

static double mc_ref_prob(const mc_aux_t *ma, double *_PD, double *f_exp)
{
	int k, i;
	long double PD = 0., Pref = 0., Ef = 0.;
	for (k = 0; k <= ma->M; ++k) {
		long double x = 1., y = 0.;
		double *bk = ma->beta + k * 3;
		for (i = 0; i < ma->n; ++i) {
			double *pdg = ma->pdg + i * 3;
			double z = pdg[0] * bk[0] + pdg[1] * bk[1] + pdg[2] * bk[2];
			x *= z;
			y += (pdg[1] * bk[1] + 2. * pdg[2] * bk[2]) / z;
		}
		PD += x * ma->alpha[k];
		Ef += x * y * ma->alpha[k];
	}
	for (k = 0; k <= ma->n * 2; ++k) {
		long double x = 1.0;
		for (i = 0; i < ma->n; ++i)
			x *= ma->pdg[i * 3 + 2] * ma->beta[k * 3 + 2];
		Pref += x * ma->alpha[k];
	}
	*f_exp = (double)(Ef / PD / ma->M);
	*_PD = PD;
	return Pref / PD;
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
// calculate z_{nr}^{(k)}
static void mc_cal_z(mc_aux_t *ma, int k)
{
	double *z[2], *tmp, *bk, *pdg;
	int i, j;
	z[0] = ma->z;
	z[1] = ma->zswap;
	bk = ma->beta + k * 3; pdg = ma->pdg;
	z[0][0] = 1.; z[0][1] = z[0][2] = 0.;
	for (j = 0; j < ma->n; ++j) {
		int max = (j + 1) * 2;
		double p[3];
		pdg = ma->pdg + j * 3;
		p[0] = bk[0] * pdg[0]; p[1] = bk[1] * pdg[1]; p[2] = bk[2] * pdg[2];
		z[1][0] = p[0] * z[0][0];
		z[1][1] = p[0] * z[0][1] + p[1] * z[0][0];
		for (i = 2; i <= max; ++i)
			z[1][i] = p[0] * z[0][i] + p[1] * z[0][i-1] + p[2] * z[0][i-2];
		if (j < ma->n - 1) z[1][max+1] = z[1][max+2] = 0.;
		tmp = z[0]; z[0] = z[1]; z[1] = tmp;
	}
	if (z[0] != ma->z) memcpy(ma->z, z[0], sizeof(double) * (2 * ma->n + 1));
}
// Warning: this is cubic in ma->n, very sloooooow
static void mc_add_afs(mc_aux_t *ma, double PD, double *f_map, double *p_map)
{
	int k, l;
	double sum = 0.;
	memset(ma->afs1, 0, sizeof(double) * (2 * ma->n + 1));
	*f_map = *p_map = -1.;
	for (k = 0; k <= 2 * ma->n; ++k) {
		mc_cal_z(ma, k);
		for (l = 0; l <= 2 * ma->n; ++l)
			ma->afs1[l] += ma->alpha[k] * ma->z[l] / PD;
	}
	for (k = 0; k <= ma->M; ++k)
		if (isnan(ma->afs1[k]) || isinf(ma->afs1[k])) return;
	for (k = 0; k <= 2 * ma->n; ++k) {
		ma->afs[k] += ma->afs1[k];
		sum += ma->afs1[k];
	}
	{
		int max_k = 0;
		double max = -1., e = 0.;
		for (k = 0; k <= 2 * ma->n; ++k) {
			if (ma->afs1[k] > max) max = ma->afs1[k], max_k = k;
			e += k * ma->afs1[k];
		}
		*f_map = .5 * max_k / ma->n; *p_map = max; // e should equal mc_rst_t::f_exp
//		printf(" * %.3lg:%.3lg:%.3lg:%.3lg * ", sum, 1.-.5*max_k/ma->n, max, 1.-.5*e/ma->n);
	}
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
	if (level >= 2) // quadratic-time calculations; necessary for genotyping
		rst->p_ref = mc_ref_prob(ma, &rst->PD, &rst->f_exp);
	if (level >= 3)
		mc_add_afs(ma, rst->PD, &rst->f_map, &rst->p_map);
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
