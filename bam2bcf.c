#include <math.h>
#include <stdint.h>
#include "bam.h"
#include "kstring.h"
#include "bam2bcf.h"
#include "errmod.h"
#include "bcftools/bcf.h"

extern	void ks_introsort_uint32_t(size_t n, uint32_t a[]);

#define CALL_ETA 0.03f
#define CALL_MAX 256
#define CALL_DEFTHETA 0.83f

#define CAP_DIST 25

struct __bcf_callaux_t {
	int max_bases, capQ, min_baseQ;
	uint16_t *bases;
	errmod_t *e;
};

bcf_callaux_t *bcf_call_init(double theta, int min_baseQ)
{
	bcf_callaux_t *bca;
	if (theta <= 0.) theta = CALL_DEFTHETA;
	bca = calloc(1, sizeof(bcf_callaux_t));
	bca->capQ = 60;
	bca->min_baseQ = min_baseQ;
	bca->e = errmod_init(1. - theta);
	return bca;
}

void bcf_call_destroy(bcf_callaux_t *bca)
{
	if (bca == 0) return;
	errmod_destroy(bca->e);
	free(bca->bases); free(bca);
}

int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base /*4-bit*/, bcf_callaux_t *bca, bcf_callret1_t *r)
{
	int i, n, ref4;
	memset(r, 0, sizeof(bcf_callret1_t));
	ref4 = bam_nt16_nt4_table[ref_base];
	if (_n == 0) return -1;

	// enlarge the bases array if necessary
	if (bca->max_bases < _n) {
		bca->max_bases = _n;
		kroundup32(bca->max_bases);
		bca->bases = (uint16_t*)realloc(bca->bases, 2 * bca->max_bases);
	}
	// fill the bases array
	memset(r, 0, sizeof(bcf_callret1_t));
	for (i = n = 0; i < _n; ++i) {
		const bam_pileup1_t *p = pl + i;
		int q, b, mapQ, baseQ, is_diff, min_dist;
		// set base
		if (p->is_del || (p->b->core.flag&BAM_FUNMAP)) continue; // skip unmapped reads and deleted bases
		baseQ = q = (int)bam1_qual(p->b)[p->qpos]; // base quality
		if (q < bca->min_baseQ) continue;
		mapQ = p->b->core.qual < bca->capQ? p->b->core.qual : bca->capQ;
		if (q > mapQ) q = mapQ;
		if (q > 63) q = 63;
		if (q < 4) q = 4;
		b = bam1_seqi(bam1_seq(p->b), p->qpos); // base
		b = bam_nt16_nt4_table[b? b : ref_base]; // b is the 2-bit base
		bca->bases[n++] = q<<5 | (int)bam1_strand(p->b)<<4 | b;
		// collect annotations
		r->qsum[b] += q;
		is_diff = (ref4 < 4 && b == ref4)? 0 : 1;
		++r->anno[0<<2|is_diff<<1|bam1_strand(p->b)];
		min_dist = p->b->core.l_qseq - 1 - p->qpos;
		if (min_dist > p->qpos) min_dist = p->qpos;
		if (min_dist > CAP_DIST) min_dist = CAP_DIST;
		r->anno[1<<2|is_diff<<1|0] += baseQ;
		r->anno[1<<2|is_diff<<1|1] += baseQ * baseQ;
		r->anno[2<<2|is_diff<<1|0] += mapQ;
		r->anno[2<<2|is_diff<<1|1] += mapQ * mapQ;
		r->anno[3<<2|is_diff<<1|0] += min_dist;
		r->anno[3<<2|is_diff<<1|1] += min_dist * min_dist;
	}
	r->depth = n;
	// glfgen
	errmod_cal(bca->e, n, 5, bca->bases, r->p);
	return r->depth;
}

int bcf_call_combine(int n, const bcf_callret1_t *calls, int ref_base /*4-bit*/, bcf_call_t *call)
{
	int ref4, i, j, qsum[4];
	int64_t tmp;
	call->ori_ref = ref4 = bam_nt16_nt4_table[ref_base];
	if (ref4 > 4) ref4 = 4;
	// calculate qsum
	memset(qsum, 0, 4 * sizeof(int));
	for (i = 0; i < n; ++i)
		for (j = 0; j < 4; ++j)
			qsum[j] += calls[i].qsum[j];
	for (j = 0; j < 4; ++j) qsum[j] = qsum[j] << 2 | j;
	// find the top 2 alleles
	for (i = 1; i < 4; ++i) // insertion sort
		for (j = i; j > 0 && qsum[j] < qsum[j-1]; --j)
			tmp = qsum[j], qsum[j] = qsum[j-1], qsum[j-1] = tmp;
	// set the reference allele and alternative allele(s)
	for (i = 0; i < 5; ++i) call->a[i] = -1;
	call->unseen = -1;
	call->a[0] = ref4;
	for (i = 3, j = 1; i >= 0; --i) {
		if ((qsum[i]&3) != ref4) {
			if (qsum[i]>>2 != 0) call->a[j++] = qsum[i]&3;
			else break;
		}
	}
	if (((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0)
		call->unseen = j, call->a[j++] = qsum[i]&3;
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
	// combine annotations
	memset(call->anno, 0, 16 * sizeof(int));
	for (i = call->depth = 0, tmp = 0; i < n; ++i) {
		call->depth += calls[i].depth;
		for (j = 0; j < 16; ++j) call->anno[j] += calls[i].anno[j];
	}
	return 0;
}

int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b)
{
	kstring_t s;
	int i;
	b->n_smpl = bc->n;
	b->tid = tid; b->pos = pos; b->qual = 0;
	s.s = b->str; s.m = b->m_str; s.l = 0;
	kputc('\0', &s);
	kputc("ACGTN"[bc->ori_ref], &s); kputc('\0', &s);
	for (i = 1; i < 5; ++i) {
		if (bc->a[i] < 0) break;
		if (i > 1) kputc(',', &s);
		kputc(bc->unseen == i? 'X' : "ACGT"[bc->a[i]], &s);
	}
	kputc('\0', &s);
	kputc('\0', &s);
	// INFO
	kputs("I16=", &s);
	for (i = 0; i < 16; ++i) {
		if (i) kputc(',', &s);
		kputw(bc->anno[i], &s);
	}
	kputc('\0', &s);
	// FMT
	kputs("PL", &s); kputc('\0', &s);
	b->m_str = s.m; b->str = s.s; b->l_str = s.l;
	bcf_sync(b);
	memcpy(b->gi[0].data, bc->PL, b->gi[0].len * bc->n);
	return 0;
}
