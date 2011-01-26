#include <math.h>
#include "errmod.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint16_t)

typedef struct __errmod_coef_t {
	double *fk, *beta, *lhet;
} errmod_coef_t;

typedef struct {
	double fsum[16], bsum[16];
	uint32_t c[16];
} call_aux_t;

static errmod_coef_t *cal_coef(double depcorr, double eta)
{
	int k, n, q;
	long double sum, sum1;
	double *lC;
	errmod_coef_t *ec;

	ec = calloc(1, sizeof(errmod_coef_t));
	// initialize ->fk
	ec->fk = (double*)calloc(256, sizeof(double));
	ec->fk[0] = 1.0;
	for (n = 1; n != 256; ++n)
		ec->fk[n] = pow(1. - depcorr, n) * (1.0 - eta) + eta;
	// initialize ->coef
	ec->beta = (double*)calloc(256 * 256 * 64, sizeof(double));
	lC = (double*)calloc(256 * 256, sizeof(double));
	for (n = 1; n != 256; ++n) {
		double lgn = lgamma(n+1);
		for (k = 1; k <= n; ++k)
			lC[n<<8|k] = lgn - lgamma(k+1) - lgamma(n-k+1);
	}
	for (q = 1; q != 64; ++q) {
		double e = pow(10.0, -q/10.0);
		double le = log(e);
		double le1 = log(1.0 - e);
		for (n = 1; n <= 255; ++n) {
			double *beta = ec->beta + (q<<16|n<<8);
			sum1 = sum = 0.0;
			for (k = n; k >= 0; --k, sum1 = sum) {
				sum = sum1 + expl(lC[n<<8|k] + k*le + (n-k)*le1);
				beta[k] = -10. / M_LN10 * logl(sum1 / sum);
			}
		}
	}
	// initialize ->lhet
	ec->lhet = (double*)calloc(256 * 256, sizeof(double));
	for (n = 0; n < 256; ++n)
		for (k = 0; k < 256; ++k)
			ec->lhet[n<<8|k] = lC[n<<8|k] - M_LN2 * n;
	free(lC);
	return ec;
}

errmod_t *errmod_init(float depcorr)
{
	errmod_t *em;
	em = (errmod_t*)calloc(1, sizeof(errmod_t));
	em->depcorr = depcorr;
	em->coef = cal_coef(depcorr, 0.03);
	return em;
}

void errmod_destroy(errmod_t *em)
{
	if (em == 0) return;
	free(em->coef->lhet); free(em->coef->fk); free(em->coef->beta);
	free(em->coef); free(em);
}
// qual:6, strand:1, base:4
int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q)
{
	call_aux_t aux;
	int i, j, k, w[32];

	if (m > m) return -1;
	memset(q, 0, m * m * sizeof(float));
	if (n == 0) return 0;
	// calculate aux.esum and aux.fsum
	if (n > 255) { // then sample 255 bases
		ks_shuffle(uint16_t, n, bases);
		n = 255;
	}
	ks_introsort(uint16_t, n, bases);
	memset(w, 0, 32 * sizeof(int));
	memset(&aux, 0, sizeof(call_aux_t));
	for (j = n - 1; j >= 0; --j) { // calculate esum and fsum
		uint16_t b = bases[j];
		int q = b>>5 < 4? 4 : b>>5;
		if (q > 63) q = 63;
		k = b&0x1f;
		aux.fsum[k&0xf] += em->coef->fk[w[k]];
		aux.bsum[k&0xf] += em->coef->fk[w[k]] * em->coef->beta[q<<16|n<<8|aux.c[k&0xf]];
		++aux.c[k&0xf];
		++w[k];
	}
	// generate likelihood
	for (j = 0; j != m; ++j) {
		float tmp1, tmp3;
		int tmp2, bar_e;
		// homozygous
		for (k = 0, tmp1 = tmp3 = 0.0, tmp2 = 0; k != m; ++k) {
			if (k == j) continue;
			tmp1 += aux.bsum[k]; tmp2 += aux.c[k]; tmp3 += aux.fsum[k];
		}
		if (tmp2) {
			bar_e = (int)(tmp1 / tmp3 + 0.499);
			if (bar_e > 63) bar_e = 63;
			q[j*m+j] = tmp1;
		}
		// heterozygous
		for (k = j + 1; k < m; ++k) {
			int cjk = aux.c[j] + aux.c[k];
			for (i = 0, tmp2 = 0, tmp1 = tmp3 = 0.0; i < m; ++i) {
				if (i == j || i == k) continue;
				tmp1 += aux.bsum[i]; tmp2 += aux.c[i]; tmp3 += aux.fsum[i];
			}
			if (tmp2) {
				bar_e = (int)(tmp1 / tmp3 + 0.499);
				if (bar_e > 63) bar_e = 63;
				q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk<<8|aux.c[k]] + tmp1;
			} else q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk<<8|aux.c[k]]; // all the bases are either j or k
		}
		for (k = 0; k != m; ++k) if (q[j*m+k] < 0.0) q[j*m+k] = 0.0;
	}
	return 0;
}
