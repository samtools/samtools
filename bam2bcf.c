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
#define DEF_MAPQ 20

#define CAP_DIST 25

bcf_callaux_t *bcf_call_init(double theta, int min_baseQ)
{
	bcf_callaux_t *bca;
	if (theta <= 0.) theta = CALL_DEFTHETA;
	bca = calloc(1, sizeof(bcf_callaux_t));
	bca->capQ = 60;
	bca->openQ = 40; bca->extQ = 20; bca->tandemQ = 100;
	bca->min_baseQ = min_baseQ;
	bca->e = errmod_init(1. - theta);
	bca->min_frac = 0.002;
	bca->min_support = 1;
	return bca;
}

void bcf_call_destroy(bcf_callaux_t *bca)
{
	if (bca == 0) return;
	errmod_destroy(bca->e);
	free(bca->bases); free(bca->inscns); free(bca);
}
/* ref_base is the 4-bit representation of the reference base. It is
 * negative if we are looking at an indel. */
int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r)
{
    static int *var_pos = NULL, nvar_pos = 0;
	int i, n, ref4, is_indel, ori_depth = 0;
	memset(r, 0, sizeof(bcf_callret1_t));
	if (ref_base >= 0) {
		ref4 = bam_nt16_nt4_table[ref_base];
		is_indel = 0;
	} else ref4 = 4, is_indel = 1;
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
		int q, b, mapQ, baseQ, is_diff, min_dist, seqQ;
		// set base
		if (p->is_del || p->is_refskip || (p->b->core.flag&BAM_FUNMAP)) continue;
		++ori_depth;
		baseQ = q = is_indel? p->aux&0xff : (int)bam1_qual(p->b)[p->qpos]; // base/indel quality
		seqQ = is_indel? (p->aux>>8&0xff) : 99;
		if (q < bca->min_baseQ) continue;
		if (q > seqQ) q = seqQ;
		mapQ = p->b->core.qual < 255? p->b->core.qual : DEF_MAPQ; // special case for mapQ==255
		mapQ = mapQ < bca->capQ? mapQ : bca->capQ;
		if (q > mapQ) q = mapQ;
		if (q > 63) q = 63;
		if (q < 4) q = 4;
		if (!is_indel) {
			b = bam1_seqi(bam1_seq(p->b), p->qpos); // base
			b = bam_nt16_nt4_table[b? b : ref_base]; // b is the 2-bit base
			is_diff = (ref4 < 4 && b == ref4)? 0 : 1;
		} else {
			b = p->aux>>16&0x3f;
			is_diff = (b != 0);
		}
		bca->bases[n++] = q<<5 | (int)bam1_strand(p->b)<<4 | b;
		// collect annotations
		if (b < 4) r->qsum[b] += q;
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
	r->depth = n; r->ori_depth = ori_depth;
	// glfgen
	errmod_cal(bca->e, n, 5, bca->bases, r->p);

    // Calculate the Variant Distance Bias (make it optional?)
    if ( nvar_pos < _n ) {
        nvar_pos = _n;
        var_pos = realloc(var_pos,sizeof(int)*nvar_pos);
    }
    int alt_dp=0, read_len=0;
    for (i=0; i<_n; i++) {
        const bam_pileup1_t *p = pl + i;
        if ( bam1_seqi(bam1_seq(p->b),p->qpos) == ref_base ) 
            continue;

        var_pos[alt_dp] = p->qpos;
        if ( (bam1_cigar(p->b)[0]&BAM_CIGAR_MASK)==4 )
            var_pos[alt_dp] -= bam1_cigar(p->b)[0]>>BAM_CIGAR_SHIFT;

        alt_dp++;
        read_len += p->b->core.l_qseq;
    }
    float mvd=0;
    int j;
    n=0;
    for (i=0; i<alt_dp; i++) {
        for (j=0; j<i; j++) {
            mvd += abs(var_pos[i] - var_pos[j]);
            n++;
        }
    }
    r->mvd[0] = n ? mvd/n : 0;
    r->mvd[1] = alt_dp;
    r->mvd[2] = alt_dp ? read_len/alt_dp : 0;

	return r->depth;
}


void calc_vdb(int n, const bcf_callret1_t *calls, bcf_call_t *call)
{
    // Variant distance bias. Samples merged by means of DP-weighted average.

    float weight=0, tot_prob=0;

    int i;
    for (i=0; i<n; i++)
    {
        int mvd      = calls[i].mvd[0];
        int dp       = calls[i].mvd[1];
        int read_len = calls[i].mvd[2];

        if ( dp<2 ) continue;

        float prob = 0;
        if ( dp==2 )
        {
            // Exact formula
            prob = (mvd==0) ? 1.0/read_len : (read_len-mvd)*2.0/read_len/read_len;
        }
        else if ( dp==3 )
        {
            // Sin, quite accurate approximation
            float mu = read_len/2.9;
            prob = mvd>2*mu ? 0 : sin(mvd*3.14/2/mu) / (4*mu/3.14);
        }
        else
        {
            // Scaled gaussian curve, crude approximation, but behaves well. Using fixed depth for bigger depths.
            if ( dp>5 )
                dp = 5;
            float sigma2 = (read_len/1.9/(dp+1)) * (read_len/1.9/(dp+1));
            float norm   = 1.125*sqrt(2*3.14*sigma2);
            float mu     = read_len/2.9;
            if ( mvd < mu )
                prob = exp(-(mvd-mu)*(mvd-mu)/2/sigma2)/norm;
            else
                prob = exp(-(mvd-mu)*(mvd-mu)/3.125/sigma2)/norm;
        }

        //fprintf(stderr,"dp=%d mvd=%d read_len=%d -> prob=%f\n", dp,mvd,read_len,prob);
        tot_prob += prob*dp;
        weight += dp;
    }
    tot_prob = weight ? tot_prob/weight : 1; 
    //fprintf(stderr,"prob=%f\n", tot_prob);
    call->vdb = tot_prob;
}

int bcf_call_combine(int n, const bcf_callret1_t *calls, int ref_base /*4-bit*/, bcf_call_t *call)
{
	int ref4, i, j, qsum[4];
	int64_t tmp;
	if (ref_base >= 0) {
		call->ori_ref = ref4 = bam_nt16_nt4_table[ref_base];
		if (ref4 > 4) ref4 = 4;
	} else call->ori_ref = -1, ref4 = 0;
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
	if (ref_base >= 0) { // for SNPs, find the "unseen" base
		if (((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0)
			call->unseen = j, call->a[j++] = qsum[i]&3;
		call->n_alleles = j;
	} else {
		call->n_alleles = j;
		if (call->n_alleles == 1) return -1; // no reliable supporting read. stop doing anything
	}
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
			for (j = 0; j <= i; ++j)
				g[z++] = call->a[j] * 5 + call->a[i];
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
//		if (ref_base < 0) fprintf(stderr, "%d,%d,%f,%d\n", call->n_alleles, x, sum_min, call->unseen);
		call->shift = (int)(sum_min + .499);
	}
	// combine annotations
	memset(call->anno, 0, 16 * sizeof(int));
	for (i = call->depth = call->ori_depth = 0, tmp = 0; i < n; ++i) {
		call->depth += calls[i].depth;
		call->ori_depth += calls[i].ori_depth;
		for (j = 0; j < 16; ++j) call->anno[j] += calls[i].anno[j];
	}

    calc_vdb(n, calls, call);

	return 0;
}

int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b, bcf_callret1_t *bcr, int is_SP,
				 const bcf_callaux_t *bca, const char *ref)
{
	extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);
	kstring_t s;
	int i, j;
	b->n_smpl = bc->n;
	b->tid = tid; b->pos = pos; b->qual = 0;
	s.s = b->str; s.m = b->m_str; s.l = 0;
	kputc('\0', &s);
	if (bc->ori_ref < 0) { // an indel
		// write REF
		kputc(ref[pos], &s);
		for (j = 0; j < bca->indelreg; ++j) kputc(ref[pos+1+j], &s);
		kputc('\0', &s);
		// write ALT
		kputc(ref[pos], &s);
		for (i = 1; i < 4; ++i) {
			if (bc->a[i] < 0) break;
			if (i > 1) {
				kputc(',', &s); kputc(ref[pos], &s);
			}
			if (bca->indel_types[bc->a[i]] < 0) { // deletion
				for (j = -bca->indel_types[bc->a[i]]; j < bca->indelreg; ++j)
					kputc(ref[pos+1+j], &s);
			} else { // insertion; cannot be a reference unless a bug
				char *inscns = &bca->inscns[bc->a[i] * bca->maxins];
				for (j = 0; j < bca->indel_types[bc->a[i]]; ++j)
					kputc("ACGTN"[(int)inscns[j]], &s);
				for (j = 0; j < bca->indelreg; ++j) kputc(ref[pos+1+j], &s);
			}
		}
		kputc('\0', &s);
	} else { // a SNP
		kputc("ACGTN"[bc->ori_ref], &s); kputc('\0', &s);
		for (i = 1; i < 5; ++i) {
			if (bc->a[i] < 0) break;
			if (i > 1) kputc(',', &s);
			kputc(bc->unseen == i? 'X' : "ACGT"[bc->a[i]], &s);
		}
		kputc('\0', &s);
	}
	kputc('\0', &s);
	// INFO
	if (bc->ori_ref < 0) kputs("INDEL;", &s);
	kputs("DP=", &s); kputw(bc->ori_depth, &s); kputs(";I16=", &s);
	for (i = 0; i < 16; ++i) {
		if (i) kputc(',', &s);
		kputw(bc->anno[i], &s);
	}
    if ( bc->vdb!=1 )
    {
        ksprintf(&s, ";VDB=%.4f", bc->vdb);
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
		int32_t *sp = is_SP? b->gi[2].data : 0;
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
