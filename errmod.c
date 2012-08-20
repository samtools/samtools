#include <math.h>
#include "errmod.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint16_t)

/* table of constants generated for given depcorr and eta */
typedef struct __errmod_coef_t {
	double *fk, *beta, *lhet;
} errmod_coef_t;

typedef struct {
	double fsum[16], bsum[16];
	uint32_t c[16];
} call_aux_t;

/* \Gamma(n) = (n-1)! */
#define lfact(n) lgamma(n+1)

/* generates a success * trials table of bionomial probability densities (log transformed) */
static double* logbinomial_table( const int n_size )
{
    /* prob distribution for binom var is p(k) = {n! \over k! (n-k)! } p^k (1-p)^{n-k} */
    /* this calcs p(k) = {log(n!) - log(k!) - log((n-k)!) */
    int k, n;
    double *logbinom = (double*)calloc(n_size * n_size, sizeof(double));
	for (n = 1; n < n_size; ++n) {
		double lfn = lfact(n);
		for (k = 1; k <= n; ++k)
			logbinom[n<<8|k] = lfn - lfact(k) - lfact(n-k);
	}
    return logbinom;
}

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

	lC = logbinomial_table( 256 );

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

/**
 * Create errmod_t object with obj.depcorr set to depcorr and initialise
 */
errmod_t *errmod_init(double depcorr)
{
	errmod_t *em;
	em = (errmod_t*)calloc(1, sizeof(errmod_t));
	em->depcorr = depcorr;
	em->coef = cal_coef(depcorr, 0.03);
	return em;
}

/**
 * Deallocate an errmod_t object
 */
void errmod_destroy(errmod_t *em)
{
	if (em == 0) return;
	free(em->coef->lhet); free(em->coef->fk); free(em->coef->beta);
	free(em->coef); free(em);
}

//
// em: error model to fit to data
// m: number of alleles across all samples
// n: number of bases observed in sample
// bases[i]: bases observed in pileup [6 bit quality|1 bit strand|4 bit base]
// q[i*m+j]: (Output) phred-scaled likelihood of each genotype (i,j)
int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q)
{
	// Aux 
	// aux.c is total count of each base observed (ignoring strand)
	call_aux_t aux;
	// Loop variables
	int i, j, k;
	// The total count of each base observed per strand
	int w[32];

	if (m > m) return -1; // FIXME: This looks like a typo?
        /* zero out q */
	memset(q, 0, m * m * sizeof(float));
	if (n == 0) return 0;
	// calculate aux.esum and aux.fsum
	if (n > 255) { // then sample 255 bases
		ks_shuffle(uint16_t, n, bases);
		n = 255;
	}
	ks_introsort(uint16_t, n, bases);
	/* zero out w and aux */
	memset(w, 0, 32 * sizeof(int));
	memset(&aux, 0, sizeof(call_aux_t));

	for (j = n - 1; j >= 0; --j) { // calculate esum and fsum
		uint16_t b = bases[j];
		/* extract quality and cap at 63 */
		int qual = b>>5 < 4? 4 : b>>5;
		if (qual > 63) qual = 63;
		/* extract base ORed with strand */ 
		int basestrand = b&0x1f;
		/* extract base */
		int base = b&0xf;
		aux.fsum[base] += em->coef->fk[w[basestrand]];
		aux.bsum[base] += em->coef->fk[w[basestrand]] * em->coef->beta[qual<<16|n<<8|aux.c[base]];
		++aux.c[base];
		++w[basestrand];
	}

	// generate likelihood
	for (j = 0; j < m; ++j) {
		float tmp1, tmp3;
		int tmp2;
		// homozygous
		for (k = 0, tmp1 = tmp3 = 0.0, tmp2 = 0; k < m; ++k) {
			if (k == j) continue;
			tmp1 += aux.bsum[k]; tmp2 += aux.c[k]; tmp3 += aux.fsum[k];
		}
		if (tmp2) {
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
				q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk<<8|aux.c[k]] + tmp1;
			} else q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk<<8|aux.c[k]]; // all the bases are either j or k
		}
		/* clamp to greater than 0 */
		for (k = 0; k < m; ++k) if (q[j*m+k] < 0.0) q[j*m+k] = 0.0;
	}

	return 0;
}
