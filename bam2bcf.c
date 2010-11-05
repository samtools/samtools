#include <math.h>
#include <stdint.h>
#include "bam.h"
#include "kstring.h"
#include "bam2bcf.h"
#include "bcftools/bcf.h"

extern	void ks_introsort_uint32_t(size_t n, uint32_t a[]);

#define CALL_ETA 0.03f
#define CALL_MAX 256
#define CALL_DEFTHETA 0.83f

#define CAP_DIST 25

bcf_callaux_t *bcf_call_init(double theta, int min_baseQ)
{
	bcf_callaux_t *bca;
	if (theta <= 0.) theta = CALL_DEFTHETA;
	bca = calloc(1, sizeof(bcf_callaux_t));
	bca->capQ = 60;
	bca->openQ = 40;
	bca->extQ = 20;
	bca->tandemQ = 100;
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
		if (p->is_del || p->is_refskip || (p->b->core.flag&BAM_FUNMAP)) continue;
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

int bcf_call_glfgen_gap(int pos, int _n, const bam_pileup1_t *pl, bcf_callaux_t *bca, bcf_callret1_t *r)
{
	int i, n, n_ins, n_del;
	memset(r, 0, sizeof(bcf_callret1_t));
	if (_n == 0) return -1;

	// enlarge the bases array if necessary
	if (bca->max_bases < _n) {
		bca->max_bases = _n;
		kroundup32(bca->max_bases);
		bca->bases = (uint16_t*)realloc(bca->bases, 2 * bca->max_bases);
	}
	// fill the bases array
	memset(r, 0, sizeof(bcf_callret1_t));
	r->indelreg = 10000;
	n_ins = n_del = 0;
	for (i = n = 0; i < _n; ++i) {
		const bam_pileup1_t *p = pl + i;
		int q, b, mapQ, indelQ, is_diff, min_dist;
		if (p->is_del || p->is_refskip || (p->b->core.flag&BAM_FUNMAP)) continue;
		{ // compute indel (base) quality
			// this can be made more efficient, but realignment is the bottleneck anyway
			int j, k, x, y, op, len = 0, max_left, max_rght, seqQ, indelreg;
			bam1_core_t *c = &p->b->core;
			uint32_t *cigar = bam1_cigar(p->b);
			uint8_t *qual = bam1_qual(p->b);
			for (k = y = 0, x = c->pos; k < c->n_cigar && y <= p->qpos; ++k) {
				op = cigar[k]&0xf;
				len = cigar[k]>>4;
				if (op == BAM_CMATCH) {
					if (pos > x && pos < x + len) break;
					x += len; y += len;
				} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += len;
				else if (op == BAM_CDEL || op == BAM_CREF_SKIP) x += len;
			}
			if (k == c->n_cigar) continue; // this actually should not happen
			max_left = max_rght = 0; indelreg = 0;
			if (pos == x + len - 1 && k+2 < c->n_cigar && ((cigar[k+1]&0xf) == BAM_CINS || (cigar[k+1]&0xf) == BAM_CDEL)
				&& (cigar[k+2]&0xf) == BAM_CMATCH)
			{
				for (j = y; j < y + len; ++j)
					if (max_left < qual[j]) max_left = qual[j];
				if ((cigar[k+1]&0xf) == BAM_CINS) y += cigar[k+1]>>4;
				else x += cigar[k+1]>>4;
				op = cigar[k+2]&0xf; len = cigar[k+2]>>4;
				for (j = y; j < y + len; ++j) {
					if (max_rght < qual[j]) max_rght = qual[j];
					if (qual[j] > BAM2BCF_INDELREG_THRES && indelreg == 0)
						indelreg = j - y + 1;
				}
				if (r->indelreg > indelreg) r->indelreg = indelreg;
			} else {
				for (j = y; j <= p->qpos; ++j)
					if (max_left < qual[j]) max_left = qual[j];
				for (j = p->qpos + 1; j < y + len; ++j)
					if (max_rght < qual[j]) max_rght = qual[j];
					
			}
			indelQ = max_left < max_rght? max_left : max_rght;
			// estimate the sequencing error rate
			seqQ = bca->openQ;
			if (p->indel != 0) seqQ += bca->extQ * (abs(p->indel) - 1); // FIXME: better to model homopolymer
			if (p->indel > 0) {
				++n_ins; r->ins_len += p->indel;
			} else if (p->indel < 0) {
				++n_del; r->del_len += -p->indel;
			}
			if (p->indel != 0) { // a different model for tandem repeats
				uint8_t *seq = bam1_seq(p->b);
				int tandemQ, qb = bam1_seqi(seq, p->qpos), l;
				for (j = p->qpos + 1; j < c->l_qseq; ++j)
					if (qb != bam1_seqi(seq, j)) break;
				l = j;
				for (j = (int)p->qpos - 1; j >= 0; --j)
					if (qb != bam1_seqi(seq, j)) break;
				l = l - (j + 1);
				tandemQ = (int)((double)(abs(p->indel)) / l * bca->tandemQ + .499);
				if (seqQ > tandemQ) seqQ = tandemQ;
			}
//			fprintf(stderr, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", bam1_qname(p->b), pos+1, p->indel, indelQ, seqQ, max_left, max_rght);
			if (indelQ > seqQ) indelQ = seqQ;
			q = indelQ;
		}
		if (q < bca->min_baseQ) continue;
		mapQ = p->b->core.qual < bca->capQ? p->b->core.qual : bca->capQ;
		if (q > mapQ) q = mapQ;
		if (q > 63) q = 63;
		if (q < 4) q = 4;
		b = p->indel? 1 : 0;
		bca->bases[n++] = q<<5 | (int)bam1_strand(p->b)<<4 | b;
		// collect annotations
		r->qsum[b] += q;
		is_diff = b;
		++r->anno[0<<2|is_diff<<1|bam1_strand(p->b)];
		min_dist = p->b->core.l_qseq - 1 - p->qpos;
		if (min_dist > p->qpos) min_dist = p->qpos;
		if (min_dist > CAP_DIST) min_dist = CAP_DIST;
		r->anno[1<<2|is_diff<<1|0] += indelQ;
		r->anno[1<<2|is_diff<<1|1] += indelQ * indelQ;
		r->anno[2<<2|is_diff<<1|0] += mapQ;
		r->anno[2<<2|is_diff<<1|1] += mapQ * mapQ;
		r->anno[3<<2|is_diff<<1|0] += min_dist;
		r->anno[3<<2|is_diff<<1|1] += min_dist * min_dist;
	}
	r->depth = n;
	r->ins_len = n_ins? r->ins_len / n_ins : 0;
	r->del_len = n_del? r->del_len / n_del : 0;
	if (r->indelreg >= 10000) r->indelreg = 0;
	// glfgen
	errmod_cal(bca->e, n, 2, bca->bases, r->p);
	return r->depth;	
}

int bcf_call_combine(int n, const bcf_callret1_t *calls, int ref_base /*4-bit*/, bcf_call_t *call)
{
	int ref4, i, j, qsum[4];
	int64_t tmp;
	call->ori_ref = ref4 = bam_nt16_nt4_table[ref_base];
	call->ins_len = call->del_len = 0; call->indelreg = 0;
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

int bcf_call_combine_gap(int n, const bcf_callret1_t *calls, bcf_call_t *call)
{
	int i, j, n_ins, n_del;
	// combine annotations
	call->ori_ref = 4;
	memset(call->anno, 0, 16 * sizeof(int));
	call->ins_len = call->del_len = 0; call->indelreg = 10000;
	for (i = call->depth = 0, n_ins = n_del = 0; i < n; ++i) {
		const bcf_callret1_t *r = calls + i;
		if (r->depth > 0) {
			call->depth += r->depth;
			if (r->ins_len > 0) {
				call->ins_len += r->ins_len;
				++n_ins;
			}
			if (r->del_len > 0) {
				call->del_len += r->del_len;
				++n_del;
			}
			if (r->indelreg > 0 && call->indelreg > r->indelreg)
				call->indelreg = r->indelreg;
			for (j = 0; j < 16; ++j) call->anno[j] += r->anno[j];
		}
	}
	if (call->depth == 0) return 0; // no indels
	call->ins_len = n_ins? call->ins_len / n_ins : 0;
	call->del_len = n_del? call->del_len / n_del : 0;
	//
	for (i = 0; i < 5; ++i) call->a[i] = -1;
	call->a[0] = 0; call->a[1] = 1;
	call->unseen = -1;
	call->n_alleles = 2;
	// set the PL array
	if (call->n < n) {
		call->n = n;
		call->PL = realloc(call->PL, 15 * n);
	}
	{
		int g[3];
		double sum_min = 0.;
		g[0] = 0; g[1] = 1; g[2] = 3;
		for (i = 0; i < n; ++i) {
			uint8_t *PL = call->PL + 3 * i;
			const bcf_callret1_t *r = calls + i;
			float min = 1e37;
			for (j = 0; j < 3; ++j)
				if (min > r->p[g[j]]) min = r->p[g[j]];
			sum_min += min;
			for (j = 0; j < 3; ++j) {
				int y;
				y = (int)(r->p[g[j]] - min + .499);
				if (y > 255) y = 255;
				PL[j] = y;
			}
		}
		call->shift = (int)(sum_min + .499);
	}
	return 0;
}

int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b, bcf_callret1_t *bcr, int is_SP)
{
	extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);
	kstring_t s;
	int i;
	b->n_smpl = bc->n;
	b->tid = tid; b->pos = pos; b->qual = 0;
	s.s = b->str; s.m = b->m_str; s.l = 0;
	kputc('\0', &s);
	if (bc->ins_len > 0 || bc->del_len > 0) { // an indel
		for (i = 0; i < bc->indelreg; ++i) kputc('N', &s);
		kputc('\0', &s);
		if (bc->ins_len > 0 && bc->del_len > 0) kputs("<INDEL>", &s);
		else if (bc->ins_len > 0) kputs("<INS>", &s);
		else if (bc->del_len > 0) kputs("<DEL>", &s);
	} else { // SNP
		kputc("ACGTN"[bc->ori_ref], &s); kputc('\0', &s);
		for (i = 1; i < 5; ++i) {
			if (bc->a[i] < 0) break;
			if (i > 1) kputc(',', &s);
			kputc(bc->unseen == i? 'X' : "ACGT"[bc->a[i]], &s);
		}
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
	kputs("PL", &s);
	if (bcr) {
		kputs(":DP", &s);
		if (is_SP) kputs(":SP", &s);
	}
	kputc('\0', &s);
	b->m_str = s.m; b->str = s.s; b->l_str = s.l;
	bcf_sync(b);
	memcpy(b->gi[0].data, bc->PL, b->gi[0].len * bc->n);
	if (bcr) {
		uint16_t *dp = (uint16_t*)b->gi[1].data;
		uint8_t *sp = is_SP? b->gi[2].data : 0;
		for (i = 0; i < bc->n; ++i) {
			bcf_callret1_t *p = bcr + i;
			dp[i] = p->depth < 0xffff? p->depth : 0xffff;
			if (is_SP) {
				if (p->anno[0] + p->anno[1] < 2 || p->anno[2] + p->anno[3] < 2
					|| p->anno[0] + p->anno[2] < 2 || p->anno[1] + p->anno[3] < 2)
				{
					sp[i] = 0;
				} else {
					double left, right, two;
					int x;
					kt_fisher_exact(p->anno[0], p->anno[1], p->anno[2], p->anno[3], &left, &right, &two);
					x = (int)(-4.343 * log(two) + .499);
					if (x > 255) x = 255;
					sp[i] = x;
				}
			}
		}
	}
	return 0;
}
