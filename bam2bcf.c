#include <math.h>
#include <stdint.h>
#include "bam.h"
#include "kstring.h"
#include "bam2bcf.h"
#include "bcf.h"

extern	void ks_introsort_uint32_t(size_t n, uint32_t a[]);

#define CALL_ETA 0.03f
#define CALL_MAX 256
#define CALL_DEFTHETA 0.85f

struct __bcf_callaux_t {
	int max_info, capQ;
	double *fk;
	uint32_t *info;
};

bcf_callaux_t *bcf_call_init(double theta)
{
	bcf_callaux_t *bca;
	int n;
	if (theta <= 0.) theta = CALL_DEFTHETA;
	bca = calloc(1, sizeof(bcf_callaux_t));
	bca->capQ = 60;
	bca->fk = calloc(CALL_MAX, sizeof(double));
	bca->fk[0] = 1.;
	for (n = 1; n < CALL_MAX; ++n)
		bca->fk[n] = theta >= 1.? 1. : pow(theta, n) * (1.0 - CALL_ETA) + CALL_ETA;
	bca->fk[CALL_MAX-1] = 0.;
	return bca;
}

void bcf_call_destroy(bcf_callaux_t *bca)
{
	if (bca) {
		free(bca->info); free(bca->fk); free(bca);
	}
}

typedef struct {
	float esum[5], fsum[5];
	uint32_t c[5];
	int w[8];
} auxaux_t;

/*
  The following code is nearly identical to bam_maqcns_glfgen() under
  the simplified SOAPsnp model. It does the following:

  1) Collect strand, base, quality and mapping quality information for
     each base and combine them in an integer:

	   x = min{baseQ,mapQ}<<24 | 1<<21 | strand<<18 | base<<16 | baseQ<<8 | mapQ

  2) Sort the list of integers for the next step.

  3) For each base, calculate e_b, the sum of weighted qualities. For
     each type of base on each strand, the best quality has the highest
     weight. Only the top 255 bases on each strand are used (different
     from maqcns).

  4) Rescale the total read depth to 255.

  5) Calculate Q(D|g) = -10\log_{10}P(D|g) (d_a is the allele count):

       Q(D|<aa>)=\sum_{b\not=a}e_b
	   Q(D|<aA>)=3*(d_a+d_A)+\sum_{b\not=a,b\not=A}e_b
 */
int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base /*4-bit*/, bcf_callaux_t *bca, bcf_callret1_t *r)
{
	int i, j, k, c, n;
	float *p = r->p;
	auxaux_t aux;

	memset(r, 0, sizeof(bcf_callret1_t));
	if (_n == 0) return -1;

	// enlarge the aux array if necessary
	if (bca->max_info < _n) {
		bca->max_info = _n;
		kroundup32(bca->max_info);
		bca->info = (uint32_t*)realloc(bca->info, 4 * bca->max_info);
	}
	// fill the aux array
	for (i = n = 0; i < _n; ++i) {
		const bam_pileup1_t *p = pl + i;
		uint32_t q, x = 0, qq;
		if (p->is_del || (p->b->core.flag&BAM_FUNMAP)) continue; // skip unmapped reads and deleted bases
		q = (uint32_t)bam1_qual(p->b)[p->qpos]; // base quality
		x |= (uint32_t)bam1_strand(p->b) << 18 | q << 8 | p->b->core.qual;
		if (p->b->core.qual < q) q = p->b->core.qual; // cap the overall quality at mapping quality
		x |= q << 24;
		qq = bam1_seqi(bam1_seq(p->b), p->qpos); // base
		q = bam_nt16_nt4_table[qq? qq : ref_base]; // q is the 2-bit base
		if (q < 4) x |= 1 << 21 | q << 16;
		
		bca->info[n++] = x;
	}
	ks_introsort_uint32_t(n, bca->info);
	r->depth = n;
	// generate esum and fsum
	memset(&aux, 0, sizeof(auxaux_t));
	for (j = n - 1, r->sum_Q2 = 0; j >= 0; --j) { // calculate esum and fsum
		uint32_t info = bca->info[j];
		int tmp;
		if (info>>24 < 4 && (info>>8&0x3f) != 0) info = 4<<24 | (info&0xffffff);
		k = info>>16&7;
		if (info>>24 > 0) {
			aux.esum[k&3] += bca->fk[aux.w[k]] * (info>>24);
			aux.fsum[k&3] += bca->fk[aux.w[k]];
			if (aux.w[k] + 1 < CALL_MAX) ++aux.w[k];
			++aux.c[k&3];
		}
		tmp = (int)(info&0xff) < bca->capQ? (int)(info&0xff) : bca->capQ;
		r->sum_Q2 += tmp * tmp;
	}
	memcpy(r->esum, aux.esum, 5 * sizeof(float));
	// rescale ->c[]
	for (j = c = 0; j != 4; ++j) c += aux.c[j];
	if (c > 255) {
		for (j = 0; j != 4; ++j) aux.c[j] = (int)(254.0 * aux.c[j] / c + 0.499);
		for (j = c = 0; j != 4; ++j) c += aux.c[j];
	}
	// generate likelihood
	for (j = 0; j != 5; ++j) {
		float tmp;
		// homozygous
		for (k = 0, tmp = 0.0; k != 5; ++k)
			if (j != k) tmp += aux.esum[k];
		p[j*5+j] = tmp; // anything that is not j
		// heterozygous
		for (k = j + 1; k < 5; ++k) {
			for (i = 0, tmp = 0.0; i != 5; ++i)
				if (i != j && i != k) tmp += aux.esum[i];
			p[j*5+k] = p[k*5+j] = 3.01 * (aux.c[j] + aux.c[k]) + tmp;
		}
	}
	return 0;
}

/*
  1) Find the top 2 bases (from esum[4]).

  2) If the reference base is among the top 2, consider the locus is
     potentially biallelic and set call->a[2] as -1; otherwise, the
     locus is potentially triallelic. If the reference is ambiguous,
     take the weakest call as the pseudo-reference.
 */
int bcf_call_combine(int n, const bcf_callret1_t *calls, int ref_base /*4-bit*/, bcf_call_t *call)
{
	int ref4, i, j;
	int64_t sum[5], tmp;
	call->ori_ref = ref4 = bam_nt16_nt4_table[ref_base];
	if (ref4 > 4) ref4 = 4;
	{ // calculate esum
		double esum[5];
		memset(esum, 0, sizeof(double) * 4);
		for (i = 0; i < n; ++i) {
			for (j = 0; j < 4; ++j)
				esum[j] += calls[i].esum[j];
		}
		for (j = 0; j < 4; ++j)
			sum[j] = (int)(esum[j] * 100. + .499) << 2 | j;
	}
	// find the top 2 alleles
	for (i = 1; i < 4; ++i) // insertion sort
		for (j = i; j > 0 && sum[j] < sum[j-1]; --j)
			tmp = sum[j], sum[j] = sum[j-1], sum[j-1] = tmp;
	// set the reference allele and alternative allele(s)
	for (i = 0; i < 5; ++i) call->a[i] = -1;
	call->unseen = -1;
	call->a[0] = ref4;
	for (i = 3, j = 1; i >= 0; --i) {
		if ((sum[i]&3) != ref4) {
			if (sum[i]>>2 != 0) call->a[j++] = sum[i]&3;
			else break;
		}
	}
	if (((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0)
		call->unseen = j, call->a[j++] = sum[i]&3;
	call->n_alleles = j;
	// set the PL array
	if (call->n < n) {
		call->n = n;
		call->PL = realloc(call->PL, 15 * n);
	}
	{
		int x, g[15], z;
		double sum_min = 0.;
		x = call->n_alleles * (call->n_alleles + 1) / 2;
		// get the possible genotypes
		for (i = z = 0; i < call->n_alleles; ++i)
			for (j = i; j < call->n_alleles; ++j)
				g[z++] = call->a[i] * 5 + call->a[j];
		for (i = 0; i < n; ++i) {
			uint8_t *PL = call->PL + x * i;
			const bcf_callret1_t *r = calls + i;
			float min = 1e37;
			for (j = 0; j < x; ++j)
				if (min > r->p[g[j]]) min = r->p[g[j]];
			sum_min += min;
			for (j = 0; j < x; ++j) {
				int y;
				y = (int)(r->p[g[j]] - min + .499);
				if (y > 255) y = 255;
				PL[j] = y;
			}
		}
		call->shift = (int)(sum_min + .499);
	}
	for (i = call->depth = 0, tmp = 0; i < n; ++i) {
		call->depth += calls[i].depth;
		tmp += calls[i].sum_Q2;
	}
	call->rmsQ = (int)(sqrt((double)tmp / call->depth) + .499);
	return 0;
}

int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b)
{
	kstring_t s;
	int i;
	b->tid = tid; b->pos = pos; b->qual = 0;
	s.s = b->str; s.m = b->m_str; s.l = 0;
	kputc('\0', &s);
	kputc("ACGTN"[bc->ori_ref], &s); kputc('\0', &s);
	for (i = 1; i < 5; ++i) {
		if (bc->a[i] < 0) break;
		if (i > 1) kputc(',', &s);
//		kputc(bc->unseen == i && i != 3? 'X' : "ACGT"[bc->a[i]], &s);
		kputc(bc->unseen == i? 'X' : "ACGT"[bc->a[i]], &s);
	}
	kputc('\0', &s);
	kputc('\0', &s);
	kputs("MQ=", &s); kputw(bc->rmsQ, &s); kputs(";DP=", &s); kputw(bc->depth, &s); kputc('\0', &s);
	kputs("PL", &s); kputc('\0', &s);
	b->m_str = s.m; b->str = s.s; b->l_str = s.l;
	bcf_sync(bc->n, b);
	memcpy(b->gi[0].data, bc->PL, b->gi[0].len * bc->n);
	return 0;
}
