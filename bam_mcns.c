#include <math.h>
#include <stdlib.h>
#include "bam_mcns.h"

#define MC_MIN_QUAL 13
#define MC_MAX_SUMQ 3000
#define MC_MAX_SUMQP 1e-300

struct __mc_aux_t {
	int n, N;
	int ref, alt;
	double theta;
	double *q2p, *pdg; // pdg -> P(D|g)
	double *alpha, *beta;
	int *qsum, *bcnt, rcnt[4];
};

void mc_init_prior(mc_aux_t *ma, double theta)
{
	double sum;
	int i;
	ma->theta = theta;
	for (i = 0, sum = 0.; i < 2 * ma->n; ++i)
		sum += (ma->alpha[i] = ma->theta / (2 * ma->n - i));
	ma->alpha[2 * ma->n] = 1. - sum;
}

mc_aux_t *mc_init(int n) // FIXME: assuming diploid
{
	mc_aux_t *ma;
	int i;
	ma = calloc(1, sizeof(mc_aux_t));
	ma->n = n; ma->N = 2 * n;
	ma->q2p = calloc(MC_MAX_SUMQ + 1, sizeof(double));
	ma->qsum = calloc(4 * ma->n, sizeof(int));
	ma->bcnt = calloc(4 * ma->n, sizeof(int));
	ma->pdg = calloc(3 * ma->n, sizeof(double));
	ma->alpha = calloc(2 * ma->n + 1, sizeof(double));
	ma->beta = calloc((2 * ma->n + 1) * 3, sizeof(double));
	for (i = 0; i <= MC_MAX_SUMQ; ++i)
		ma->q2p[i] = pow(10., -i / 10.);
	for (i = 0; i <= 2 * ma->n; ++i) {
		double *bi = ma->beta + 3 * i;
		double f = (double)i / (2 * ma->n);
		bi[0] = (1. - f) * (1. - f);
		bi[1] = 2 * f * (1. - f);
		bi[2] = f * f;
	}
	mc_init_prior(ma, 1e-3); // the simplest prior
	return ma;
}

void mc_destroy(mc_aux_t *ma)
{
	if (ma) {
		free(ma->qsum); free(ma->bcnt);
		free(ma->q2p); free(ma->pdg);
		free(ma->alpha); free(ma->beta);
		free(ma);
	}
}

static void sum_err(int *n, const bam_pileup1_t **plp, mc_aux_t *ma)
{
	int i, j;
	memset(ma->qsum, 0, sizeof(int) * 4 * ma->n);
	memset(ma->bcnt, 0, sizeof(int) * 4 * ma->n);
	memset(ma->rcnt, 0, sizeof(int) * 4);
	for (j = 0; j < ma->n; ++j) {
		int tmp, *qsum = ma->qsum + j * 4;
		int *bcnt = ma->bcnt + j * 4;
		double r = drand48(), rr;
		for (i = tmp = 0; i < n[j]; ++i) {
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
			++tmp;
		}
		if (tmp) {
			for (i = 0, rr = 0.; i < 4; ++i) {
				rr += (double)bcnt[i] / tmp;
				if (rr >= r) break;
			}
			if (i < 4) ++ma->rcnt[i];
		}
	}
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
	ma->ref = sum[3]&3; ma->alt = sum[2]&3;
	if (ref == ma->alt) tmp = ma->ref, ma->ref = ma->alt, ma->alt = tmp;
	// note that ma->ref might not be ref in case of triallele
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
// return the reference allele frequency
double mc_freq0(int ref, int *n, const bam_pileup1_t **plp, mc_aux_t *ma, int *_ref, int *alt)
{
	int i, acnt[4], j;
	double fr = -1., f0 = -1.;
	sum_err(n, plp, ma);
	set_allele(ref, ma);
	cal_pdg(ma);
	acnt[0] = acnt[1] = acnt[2] = acnt[3] = 0;
	for (i = 0; i < ma->n; ++i)
		for (j = 0; j < 4; ++j)
			acnt[j] += ma->bcnt[i * 4 + j];
	if (ma->rcnt[ma->ref] + ma->rcnt[ma->alt]) // never used...
		fr = (double)ma->rcnt[ref] / (ma->rcnt[ma->ref] + ma->rcnt[ma->alt]);
	if (acnt[ma->ref] + acnt[ma->alt])
		f0 = (double)acnt[ref] / (acnt[ma->ref] + acnt[ma->alt]);
	*_ref = ma->ref; *alt = ma->alt;
	return ma->n == 1 || fr < 0.? f0 : fr;
}
// f0 is the reference allele frequency
double mc_freq_iter(double f0, mc_aux_t *ma)
{
	double f, f3[3];
	int i, j;
	f3[0] = (1.-f0)*(1.-f0); f3[1] = 2.*f0*(1.-f0); f3[2] = f0*f0;
	for (i = 0, f = 0.; i < ma->n; ++i) {
		double up, dn, *pdg;
		pdg = ma->pdg + i * 3;
		for (j = 1, up = 0.; j < 3; ++j)
			up += j * pdg[j] * f3[j];
		for (j = 0, dn = 0.; j < 3; ++j)
			dn += pdg[j] * f3[j];
		f += up / dn;
	}
	f /= ma->n * 2.;
	return f;
}

double mc_ref_prob(mc_aux_t *ma)
{
	int k, i, g;
	long double PD = 0., Pref = 0.;
	for (k = 0; k <= ma->n * 2; ++k) {
		long double x = 1.;
		double *bk = ma->beta + k * 3;
		for (i = 0; i < ma->n; ++i) {
			double y = 0., *pdg = ma->pdg + i * 3;
			for (g = 0; g < 3; ++g)
				y += pdg[g] * bk[g];
			x *= y;
		}
		PD += x * ma->alpha[k];
	}
	for (k = 0; k <= ma->n * 2; ++k) {
		long double x = 1.0;
		for (i = 0; i < ma->n; ++i)
			x *= ma->pdg[i * 3 + 2] * ma->beta[k * 3 + 2];
		Pref += x * ma->alpha[k];
	}
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
