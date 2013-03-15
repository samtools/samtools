#include <math.h>
#include <stdint.h>
#include <assert.h>
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
    bca->per_sample_flt = 0;
    bca->npos = 100;
    bca->ref_pos = calloc(bca->npos, sizeof(int));
    bca->alt_pos = calloc(bca->npos, sizeof(int));
 	return bca;
}


static int get_position(const bam_pileup1_t *p, int *len)
{
    int icig, n_tot_bases = 0, iread = 0, edist = p->qpos + 1;
    for (icig=0; icig<p->b->core.n_cigar; icig++) 
    {
        // Conversion from uint32_t to MIDNSHP
        //  0123456
        //  MIDNSHP
        int cig  = bam1_cigar(p->b)[icig] & BAM_CIGAR_MASK;
        int ncig = bam1_cigar(p->b)[icig] >> BAM_CIGAR_SHIFT;
        if ( cig==0 )
        {
            n_tot_bases += ncig;
            iread += ncig;
        }
        else if ( cig==1 )
        {
            n_tot_bases += ncig;
            iread += ncig;
        }
        else if ( cig==4 )
        {
            iread += ncig;
            if ( iread<=p->qpos ) edist -= ncig;
        }
    }
    *len = n_tot_bases;
    return edist;
}

void bcf_call_destroy(bcf_callaux_t *bca)
{
	if (bca == 0) return;
	errmod_destroy(bca->e);
    if (bca->npos) { free(bca->ref_pos); free(bca->alt_pos); bca->npos = 0; }
	free(bca->bases); free(bca->inscns); free(bca);
}
/* ref_base is the 4-bit representation of the reference base. It is
 * negative if we are looking at an indel. */
int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r)
{
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
	for (i = n = r->n_supp = 0; i < _n; ++i) {
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
		if (is_diff) ++r->n_supp;
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

        // collect read positions for ReadPosBias
        int len, pos = get_position(p, &len);
        int epos = (double)pos/(len+1) * bca->npos;
        if ( bam1_seqi(bam1_seq(p->b),p->qpos) == ref_base )
            bca->ref_pos[epos]++;
        else
            bca->alt_pos[epos]++;
	}
	r->depth = n; r->ori_depth = ori_depth;
	// glfgen
	errmod_cal(bca->e, n, 5, bca->bases, r->p);
	return r->depth;
}

double mann_whitney_1947(int n, int m, int U)
{
    if (U<0) return 0;
    if (n==0||m==0) return U==0 ? 1 : 0;
    return (double)n/(n+m)*mann_whitney_1947(n-1,m,U-m) + (double)m/(n+m)*mann_whitney_1947(n,m-1,U);
}

void calc_ReadPosBias(bcf_callaux_t *bca, bcf_call_t *call)
{
    int i, nref = 0, nalt = 0;
    unsigned long int U = 0;
    for (i=0; i<bca->npos; i++) 
    {
        nref += bca->ref_pos[i];
        nalt += bca->alt_pos[i];
        U += nref*bca->alt_pos[i];
        bca->ref_pos[i] = 0;
        bca->alt_pos[i] = 0;
    }
#if 0
//todo
    double var = 0, avg = (double)(nref+nalt)/bca->npos;
    for (i=0; i<bca->npos; i++) 
    {
        double ediff = bca->ref_pos[i] + bca->alt_pos[i] - avg;
        var += ediff*ediff;
        bca->ref_pos[i] = 0;
        bca->alt_pos[i] = 0;
    }
    call->read_pos.avg = avg;
    call->read_pos.var = sqrt(var/bca->npos);
    call->read_pos.dp  = nref+nalt;
#endif
    if ( !nref || !nalt )
    {
        call->read_pos_bias = -1;
        return;
    }

    if ( nref>=8 || nalt>=8 )
    {
        // normal approximation
        double mean = ((double)nref*nalt+1.0)/2.0;
        double var2 = (double)nref*nalt*(nref+nalt+1.0)/12.0;
        double z    = (U-mean)/sqrt(var2);
        call->read_pos_bias = z;
        //fprintf(stderr,"nref=%d  nalt=%d  U=%ld  mean=%e  var=%e  zval=%e\n", nref,nalt,U,mean,sqrt(var2),call->read_pos_bias);
    }
    else
    {
        double p = mann_whitney_1947(nalt,nref,U);
        // biased form claimed by GATK to behave better empirically
        // double var2 = (1.0+1.0/(nref+nalt+1.0))*(double)nref*nalt*(nref+nalt+1.0)/12.0;
        double var2 = (double)nref*nalt*(nref+nalt+1.0)/12.0;
        double z;
        if ( p >= 1./sqrt(var2*2*M_PI) ) z = 0;   // equal to mean
        else
        {
            if ( U >= nref*nalt/2. ) z = sqrt(-2*log(sqrt(var2*2*M_PI)*p));
            else z = -sqrt(-2*log(sqrt(var2*2*M_PI)*p));
        }
        call->read_pos_bias = z;
        //fprintf(stderr,"nref=%d  nalt=%d  U=%ld  p=%e var2=%e  zval=%e\n", nref,nalt,U, p,var2,call->read_pos_bias);
    }
}

float mean_diff_to_prob(float mdiff, int dp, int readlen)
{
    if ( dp==2 )
    {
        if ( mdiff==0 )
            return (2.0*readlen + 4.0*(readlen-1.0))/((float)readlen*readlen);
        else
            return 8.0*(readlen - 4.0*mdiff)/((float)readlen*readlen);
    }

    // This is crude empirical approximation and is not very accurate for
    // shorter read lengths (<100bp). There certainly is a room for
    // improvement.
    const float mv[24][2] = { {0,0}, {0,0}, {0,0},
        { 9.108, 4.934}, { 9.999, 3.991}, {10.273, 3.485}, {10.579, 3.160},
        {10.828, 2.889}, {11.014, 2.703}, {11.028, 2.546}, {11.244, 2.391},
        {11.231, 2.320}, {11.323, 2.138}, {11.403, 2.123}, {11.394, 1.994},
        {11.451, 1.928}, {11.445, 1.862}, {11.516, 1.815}, {11.560, 1.761},
        {11.544, 1.728}, {11.605, 1.674}, {11.592, 1.652}, {11.674, 1.613},
        {11.641, 1.570} };

    float m, v;
    if ( dp>=24 )
    {
        m = readlen/8.;
        if (dp>100) dp = 100;
        v = 1.476/(0.182*pow(dp,0.514));
        v = v*(readlen/100.);
    }
    else
    {
        m = mv[dp][0];
        v = mv[dp][1];
        m = m*readlen/100.;
        v = v*readlen/100.;
        v *= 1.2;   // allow more variability
    }
    return 1.0/(v*sqrt(2*M_PI)) * exp(-0.5*((mdiff-m)/v)*((mdiff-m)/v));
}

void calc_vdb(bcf_callaux_t *bca, bcf_call_t *call)
{
    int i, dp = 0;
    float mean_pos = 0, mean_diff = 0;
    for (i=0; i<bca->npos; i++)
    {
        if ( !bca->alt_pos[i] ) continue;
        dp += bca->alt_pos[i];
        int j = i<bca->npos/2 ? i : bca->npos - i;
        mean_pos += bca->alt_pos[i]*j;
    }
    if ( dp<2 )
    {
        call->vdb = -1;
        return;
    }
    mean_pos /= dp;
    for (i=0; i<bca->npos; i++)
    {
        if ( !bca->alt_pos[i] ) continue;
        int j = i<bca->npos/2 ? i : bca->npos - i;
        mean_diff += bca->alt_pos[i] * fabs(j - mean_pos);
    }
    mean_diff /= dp;
    call->vdb = mean_diff_to_prob(mean_diff, dp, bca->npos);
}

/**
 *  bcf_call_combine() - sets the PL array and VDB, RPB annotations, finds the top two alleles
 *  @n:         number of samples
 *  @calls:     each sample's calls
 *  @bca:       auxiliary data structure for holding temporary values
 *  @ref_base:  the reference base
 *  @call:      filled with the annotations
 */
int bcf_call_combine(int n, const bcf_callret1_t *calls, bcf_callaux_t *bca, int ref_base /*4-bit*/, bcf_call_t *call)
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
    int qsum_tot=0;
    for (j=0; j<4; j++) { qsum_tot += qsum[j]; call->qsum[j] = 0; }
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
			if (qsum[i]>>2 != 0) 
            {
                if ( j<4 ) call->qsum[j] = (float)(qsum[i]>>2)/qsum_tot; // ref N can make j>=4
                call->a[j++]  = qsum[i]&3;
            }
			else break;
		}
        else 
            call->qsum[0] = (float)(qsum[i]>>2)/qsum_tot;
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

    calc_vdb(bca, call);
    calc_ReadPosBias(bca, call);

	return 0;
}

int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b, bcf_callret1_t *bcr, int fmt_flag,
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
	if (bc->ori_ref < 0) ksprintf(&s,"INDEL;IS=%d,%f;", bca->max_support, bca->max_frac);
	kputs("DP=", &s); kputw(bc->ori_depth, &s); kputs(";I16=", &s);
	for (i = 0; i < 16; ++i) {
		if (i) kputc(',', &s);
		kputw(bc->anno[i], &s);
	}
    //ksprintf(&s,";RPS=%d,%f,%f", bc->read_pos.dp,bc->read_pos.avg,bc->read_pos.var);
    ksprintf(&s,";QS=%f,%f,%f,%f", bc->qsum[0],bc->qsum[1],bc->qsum[2],bc->qsum[3]);
    if (bc->vdb != -1)
        ksprintf(&s, ";VDB=%e", bc->vdb);
    if (bc->read_pos_bias != -1 )
        ksprintf(&s, ";RPB=%e", bc->read_pos_bias);
	kputc('\0', &s);
	// FMT
	kputs("PL", &s);
	if (bcr && fmt_flag) {
		if (fmt_flag & B2B_FMT_DP) kputs(":DP", &s);
		if (fmt_flag & B2B_FMT_DV) kputs(":DV", &s);
		if (fmt_flag & B2B_FMT_SP) kputs(":SP", &s);
	}
	kputc('\0', &s);
	b->m_str = s.m; b->str = s.s; b->l_str = s.l;
	bcf_sync(b);
	memcpy(b->gi[0].data, bc->PL, b->gi[0].len * bc->n);
	if (bcr && fmt_flag) {
		uint16_t *dp = (fmt_flag & B2B_FMT_DP)? b->gi[1].data : 0;
		uint16_t *dv = (fmt_flag & B2B_FMT_DV)? b->gi[1 + ((fmt_flag & B2B_FMT_DP) != 0)].data : 0;
		int32_t  *sp = (fmt_flag & B2B_FMT_SP)? b->gi[1 + ((fmt_flag & B2B_FMT_DP) != 0) + ((fmt_flag & B2B_FMT_DV) != 0)].data : 0;
		for (i = 0; i < bc->n; ++i) {
			bcf_callret1_t *p = bcr + i;
			if (dp) dp[i] = p->depth  < 0xffff? p->depth  : 0xffff;
			if (dv) dv[i] = p->n_supp < 0xffff? p->n_supp : 0xffff;
			if (sp) {
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
