/*  bam2bcf.c -- variant calling.

    Copyright (C) 2010-2012 Broad Institute.
    Copyright (C) 2012-2014 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <float.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>
#include <htslib/kfunc.h>
#include "bam2bcf.h"

extern  void ks_introsort_uint32_t(size_t n, uint32_t a[]);

#define CALL_DEFTHETA 0.83
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
    bca->ref_pos = malloc(bca->npos*sizeof(int));
    bca->alt_pos = malloc(bca->npos*sizeof(int));
    bca->nqual = 60;
    bca->ref_mq  = malloc(bca->nqual*sizeof(int));
    bca->alt_mq  = malloc(bca->nqual*sizeof(int));
    bca->ref_bq  = malloc(bca->nqual*sizeof(int));
    bca->alt_bq  = malloc(bca->nqual*sizeof(int));
    bca->fwd_mqs = malloc(bca->nqual*sizeof(int));
    bca->rev_mqs = malloc(bca->nqual*sizeof(int));
    return bca;
}

void bcf_call_destroy(bcf_callaux_t *bca)
{
    if (bca == 0) return;
    errmod_destroy(bca->e);
    if (bca->npos) { free(bca->ref_pos); free(bca->alt_pos); bca->npos = 0; }
    free(bca->ref_mq); free(bca->alt_mq); free(bca->ref_bq); free(bca->alt_bq);
    free(bca->fwd_mqs); free(bca->rev_mqs);
    bca->nqual = 0;
    free(bca->bases); free(bca->inscns); free(bca);
}

// position in the sequence with respect to the aligned part of the read
static int get_position(const bam_pileup1_t *p, int *len)
{
    int icig, n_tot_bases = 0, iread = 0, edist = p->qpos + 1;
    for (icig=0; icig<p->b->core.n_cigar; icig++)
    {
        int cig  = bam_get_cigar(p->b)[icig] & BAM_CIGAR_MASK;
        int ncig = bam_get_cigar(p->b)[icig] >> BAM_CIGAR_SHIFT;
        if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
        {
            n_tot_bases += ncig;
            iread += ncig;
            continue;
        }
        if ( cig==BAM_CINS )
        {
            n_tot_bases += ncig;
            iread += ncig;
            continue;
        }
        if ( cig==BAM_CSOFT_CLIP )
        {
            iread += ncig;
            if ( iread<=p->qpos ) edist -= ncig;
            continue;
        }
        if ( cig==BAM_CDEL ) continue;
        if ( cig==BAM_CHARD_CLIP ) continue;
        if ( cig==BAM_CPAD ) continue;
        if ( cig==BAM_CREF_SKIP ) continue;
        fprintf(stderr,"todo: cigar %d\n", cig);
        assert(0);
    }
    *len = n_tot_bases;
    return edist;
}

void bcf_callaux_clean(bcf_callaux_t *bca, bcf_call_t *call)
{
    memset(bca->ref_pos,0,sizeof(int)*bca->npos);
    memset(bca->alt_pos,0,sizeof(int)*bca->npos);
    memset(bca->ref_mq,0,sizeof(int)*bca->nqual);
    memset(bca->alt_mq,0,sizeof(int)*bca->nqual);
    memset(bca->ref_bq,0,sizeof(int)*bca->nqual);
    memset(bca->alt_bq,0,sizeof(int)*bca->nqual);
    memset(bca->fwd_mqs,0,sizeof(int)*bca->nqual);
    memset(bca->rev_mqs,0,sizeof(int)*bca->nqual);
    if ( call->ADF ) memset(call->ADF,0,sizeof(int32_t)*(call->n+1)*B2B_MAX_ALLELES);
    if ( call->ADR ) memset(call->ADR,0,sizeof(int32_t)*(call->n+1)*B2B_MAX_ALLELES);
}

/*
    Notes:
    - Called from bam_plcmd.c by mpileup. Amongst other things, sets the bcf_callret1_t.qsum frequencies
        which are carried over via bcf_call_combine and bcf_call2bcf to the output BCF as the QS annotation.
        Later it's used for multiallelic calling by bcftools -m
    - ref_base is the 4-bit representation of the reference base. It is negative if we are looking at an indel.
 */
/*
 * This function is called once for each sample.
 * _n is number of pilesups pl contributing reads to this sample
 * pl is pointer to array of _n pileups (one pileup per read)
 * ref_base is the 4-bit representation of the reference base. It is negative if we are looking at an indel.
 * bca is the settings to perform calls across all samples
 * r is the returned value of the call
 */
int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r)
{
    int i, n, ref4, is_indel, ori_depth = 0;

    // clean from previous run
    r->ori_depth = 0;
    r->mq0 = 0;
    memset(r->qsum,0,sizeof(float)*4);
    memset(r->anno,0,sizeof(double)*16);
    memset(r->p,0,sizeof(float)*25);

    if (ref_base >= 0) {
        ref4 = seq_nt16_int[ref_base];
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
    for (i = n = 0; i < _n; ++i) {
        const bam_pileup1_t *p = pl + i;
        int q, b, mapQ, baseQ, is_diff, min_dist, seqQ;
        // set base
        if (p->is_del || p->is_refskip || (p->b->core.flag&BAM_FUNMAP)) continue;
        ++ori_depth;
        mapQ  = p->b->core.qual < 255? p->b->core.qual : DEF_MAPQ; // special case for mapQ==255
        if ( !mapQ ) r->mq0++;
        baseQ = q = is_indel? p->aux&0xff : (int)bam_get_qual(p->b)[p->qpos]; // base/indel quality
        seqQ = is_indel? (p->aux>>8&0xff) : 99;
        if (q < bca->min_baseQ) continue;
        if (q > seqQ) q = seqQ;
        mapQ = mapQ < bca->capQ? mapQ : bca->capQ;
        if (q > mapQ) q = mapQ;
        if (q > 63) q = 63;
        if (q < 4) q = 4;       // MQ=0 reads count as BQ=4
        if (!is_indel) {
            b = bam_seqi(bam_get_seq(p->b), p->qpos); // base
            b = seq_nt16_int[b? b : ref_base]; // b is the 2-bit base
            is_diff = (ref4 < 4 && b == ref4)? 0 : 1;
        } else {
            b = p->aux>>16&0x3f;
            is_diff = (b != 0);
        }
        bca->bases[n++] = q<<5 | (int)bam_is_rev(p->b)<<4 | b;
        // collect annotations
        if (b < 4)
        {
            r->qsum[b] += q;
            if ( r->ADF )
            {
                if ( bam_is_rev(p->b) )
                    r->ADR[b]++;
                else
                    r->ADF[b]++;
            }
        }
        ++r->anno[0<<2|is_diff<<1|bam_is_rev(p->b)];
        min_dist = p->b->core.l_qseq - 1 - p->qpos;
        if (min_dist > p->qpos) min_dist = p->qpos;
        if (min_dist > CAP_DIST) min_dist = CAP_DIST;
        r->anno[1<<2|is_diff<<1|0] += baseQ;
        r->anno[1<<2|is_diff<<1|1] += baseQ * baseQ;
        r->anno[2<<2|is_diff<<1|0] += mapQ;
        r->anno[2<<2|is_diff<<1|1] += mapQ * mapQ;
        r->anno[3<<2|is_diff<<1|0] += min_dist;
        r->anno[3<<2|is_diff<<1|1] += min_dist * min_dist;

        // collect for bias tests
        if ( baseQ > 59 ) baseQ = 59;
        if ( mapQ > 59 ) mapQ = 59;
        int len, pos = get_position(p, &len);
        int epos = (double)pos/(len+1) * bca->npos;
        int ibq  = baseQ/60. * bca->nqual;
        int imq  = mapQ/60. * bca->nqual;
        if ( bam_is_rev(p->b) ) bca->rev_mqs[imq]++;
        else bca->fwd_mqs[imq]++;
        if ( bam_seqi(bam_get_seq(p->b),p->qpos) == ref_base )
        {
            bca->ref_pos[epos]++;
            bca->ref_bq[ibq]++;
            bca->ref_mq[imq]++;
        }
        else
        {
            bca->alt_pos[epos]++;
            bca->alt_bq[ibq]++;
            bca->alt_mq[imq]++;
        }
    }
    r->ori_depth = ori_depth;
    // glfgen
    errmod_cal(bca->e, n, 5, bca->bases, r->p); // calculate PL of each genotype
    return n;
}


/*
 *  calc_vdb() - returns value between zero (most biased) and one (no bias)
 *               on success, or HUGE_VAL when VDB cannot be calculated because
 *               of insufficient depth (<2x)
 *
 *  Variant Distance Bias tests if the variant bases are positioned within the
 *  reads with sufficient randomness. Unlike other tests, it looks only at
 *  variant reads and therefore gives different kind of information than Read
 *  Position Bias for instance. VDB was developed for detecting artefacts in
 *  RNA-seq calls where reads from spliced transcripts span splice site
 *  boundaries.  The current implementation differs somewhat from the original
 *  version described in supplementary material of PMID:22524474, but the idea
 *  remains the same. (Here the random variable tested is the average distance
 *  from the averaged position, not the average pairwise distance.)
 *
 *  For coverage of 2x, the calculation is exact but is approximated for the
 *  rest. The result is most accurate between 4-200x. For 3x or >200x, the
 *  reported values are slightly more favourable than those of a true random
 *  distribution.
 */
double calc_vdb(int *pos, int npos)
{
    // Note well: the parameters were obtained by fitting to simulated data of
    // 100bp reads. This assumes rescaling to 100bp in bcf_call_glfgen().
    const int readlen = 100;
    assert( npos==readlen );

    #define nparam 15
    const float param[nparam][3] = { {3,0.079,18}, {4,0.09,19.8}, {5,0.1,20.5}, {6,0.11,21.5},
        {7,0.125,21.6}, {8,0.135,22}, {9,0.14,22.2}, {10,0.153,22.3}, {15,0.19,22.8},
        {20,0.22,23.2}, {30,0.26,23.4}, {40,0.29,23.5}, {50,0.35,23.65}, {100,0.5,23.7},
        {200,0.7,23.7} };

    int i, dp = 0;
    float mean_pos = 0, mean_diff = 0;
    for (i=0; i<npos; i++)
    {
        if ( !pos[i] ) continue;
        dp += pos[i];
        mean_pos += pos[i]*i;
    }
    if ( dp<2 ) return HUGE_VAL;     // one or zero reads can be placed anywhere

    mean_pos /= dp;
    for (i=0; i<npos; i++)
    {
        if ( !pos[i] ) continue;
        mean_diff += pos[i] * fabs(i - mean_pos);
    }
    mean_diff /= dp;

    int ipos = mean_diff;   // tuned for float-to-int implicit conversion
    if ( dp==2 )
        return (2*readlen-2*(ipos+1)-1)*(ipos+1)/(readlen-1)/(readlen*0.5);

    if ( dp>=200 )
        i = nparam; // shortcut for big depths
    else
    {
        for (i=0; i<nparam; i++)
            if ( param[i][0]>=dp ) break;
    }
    float pshift, pscale;
    if ( i==nparam )
    {
        // the depth is too high, go with 200x
        pscale = param[nparam-1][1];
        pshift = param[nparam-1][2];
    }
    else if ( i>0 && param[i][0]!=dp )
    {
        // linear interpolation of parameters
        pscale = (param[i-1][1] + param[i][1])*0.5;
        pshift = (param[i-1][2] + param[i][2])*0.5;
    }
    else
    {
        pscale = param[i][1];
        pshift = param[i][2];
    }
    return 0.5*kf_erfc(-(mean_diff-pshift)*pscale);
}

double calc_chisq_bias(int *a, int *b, int n)
{
    int na = 0, nb = 0, i, ndf = n;
    for (i=0; i<n; i++) na += a[i];
    for (i=0; i<n; i++) nb += b[i];
    if ( !na || !nb ) return HUGE_VAL;

    double chisq = 0;
    for (i=0; i<n; i++)
    {
        if ( !a[i] && !b[i] ) ndf--;
        else
        {
            double tmp = a[i] - b[i];
            chisq += tmp*tmp/(a[i]+b[i]);
        }
    }
    /*
        kf_gammq: incomplete gamma function Q(a,x) = 1 - P(a,x) = Gamma(a,x)/Gamma(a)
        1 if the distributions are identical, 0 if very different
    */
    double prob = kf_gammaq(0.5*ndf, 0.5*chisq);
    return prob;
}

double mann_whitney_1947(int n, int m, int U)
{
    if (U<0) return 0;
    if (n==0||m==0) return U==0 ? 1 : 0;
    return (double)n/(n+m)*mann_whitney_1947(n-1,m,U-m) + (double)m/(n+m)*mann_whitney_1947(n,m-1,U);
}

double mann_whitney_1947_cdf(int n, int m, int U)
{
    int i;
    double sum = 0;
    for (i=0; i<=U; i++)
        sum += mann_whitney_1947(n,m,i);
    return sum;
}

double calc_mwu_bias_cdf(int *a, int *b, int n)
{
    int na = 0, nb = 0, i;
    double U = 0, ties = 0;
    for (i=0; i<n; i++)
    {
        na += a[i];
        U  += a[i] * (nb + b[i]*0.5);
        nb += b[i];
        if ( a[i] && b[i] )
        {
            double tie = a[i] + b[i];
            ties += (tie*tie-1)*tie;
        }
    }
    if ( !na || !nb ) return HUGE_VAL;

    // Always work with the smaller U
    double U_min = ((double)na * nb) - U;
    if ( U < U_min ) U_min = U;

    if ( na==1 ) return 2.0 * (floor(U_min)+1) / (nb+1);
    if ( nb==1 ) return 2.0 * (floor(U_min)+1) / (na+1);

    // Normal approximation, very good for na>=8 && nb>=8 and reasonable if na<8 or nb<8
    if ( na>=8 || nb>=8 )
    {
        double mean = ((double)na*nb)*0.5;
        // Correction for ties:
        //      double N = na+nb;
        //      double var2 = (N*N-1)*N-ties;
        //      if ( var2==0 ) return 1.0;
        //      var2 *= ((double)na*nb)/N/(N-1)/12.0;
        // No correction for ties:
        double var2 = ((double)na*nb)*(na+nb+1)/12.0;
        double z = (U_min - mean)/sqrt(2*var2);   // z is N(0,1)
        return 2.0 - kf_erfc(z);  // which is 1 + erf(z)
    }

    // Exact calculation
    double pval = 2*mann_whitney_1947_cdf(na,nb,U_min);
    return pval>1 ? 1 : pval;
}

double calc_mwu_bias(int *a, int *b, int n)
{
    int na = 0, nb = 0, i;
    double U = 0, ties = 0;
    for (i=0; i<n; i++)
    {
        na += a[i];
        U  += a[i] * (nb + b[i]*0.5);
        nb += b[i];
        if ( a[i] && b[i] )
        {
            double tie = a[i] + b[i];
            ties += (tie*tie-1)*tie;
        }
    }
    if ( !na || !nb ) return HUGE_VAL;
    if ( na==1 || nb==1 ) return 1.0;       // Flat probability, all U values are equally likely

    double mean = ((double)na*nb)*0.5;
    if ( na==2 || nb==2 )
    {
        // Linear approximation
        return U>mean ? (2.0*mean-U)/mean : U/mean;
    }
    // Correction for ties:
    //      double N = na+nb;
    //      double var2 = (N*N-1)*N-ties;
    //      if ( var2==0 ) return 1.0;
    //      var2 *= ((double)na*nb)/N/(N-1)/12.0;
    // No correction for ties:
    double var2 = ((double)na*nb)*(na+nb+1)/12.0;
    if ( na>=8 || nb>=8 )
    {
        // Normal approximation, very good for na>=8 && nb>=8 and reasonable if na<8 or nb<8
        return exp(-0.5*(U-mean)*(U-mean)/var2);
    }

    // Exact calculation
    return mann_whitney_1947(na,nb,U) * sqrt(2*M_PI*var2);
}

static inline double logsumexp2(double a, double b)
{
    if ( a>b )
        return log(1 + exp(b-a)) + a;
    else
        return log(1 + exp(a-b)) + b;
}

void calc_SegBias(const bcf_callret1_t *bcr, bcf_call_t *call)
{
    call->seg_bias = HUGE_VAL;
    if ( !bcr ) return;

    int nr = call->anno[2] + call->anno[3]; // number of observed non-reference reads
    if ( !nr ) return;

    int avg_dp = (call->anno[0] + call->anno[1] + nr) / call->n;    // average depth
    double M   = floor((double)nr / avg_dp + 0.5);   // an approximate number of variants samples in the population
    if ( M>call->n ) M = call->n;       // clamp M at the number of samples
    else if ( M==0 ) M = 1;
    double f = M / 2. / call->n;        // allele frequency
    double p = (double) nr / call->n;   // number of variant reads per sample expected if variant not real (poisson)
    double q = (double) nr / M;         // number of variant reads per sample expected if variant is real (poisson)
    double sum = 0;
    const double log2 = log(2.0);

    // fprintf(stderr,"M=%.1f  p=%e q=%e f=%f  dp=%d\n",M,p,q,f,avg_dp);
    int i;
    for (i=0; i<call->n; i++)
    {
        int oi = bcr[i].anno[2] + bcr[i].anno[3];       // observed number of non-ref reads
        double tmp;
        if ( oi )
        {
            // tmp = log(f) + oi*log(q/p) - q + log(2*(1-f) + f*pow(2,oi)*exp(-q)) + p; // this can under/overflow
            tmp = logsumexp2(log(2*(1-f)), log(f) + oi*log2 - q);
            tmp += log(f) + oi*log(q/p) - q + p;
        }
        else
            tmp = log(2*f*(1-f)*exp(-q) + f*f*exp(-2*q) + (1-f)*(1-f)) + p;
        sum += tmp;
        // fprintf(stderr,"oi=%d %e\n", oi,tmp);
    }
    call->seg_bias = sum;
}

/**
 *  bcf_call_combine() - sets the PL array and VDB, RPB annotations, finds the top two alleles
 *  @n:         number of samples
 *  @calls:     each sample's calls
 *  @bca:       auxiliary data structure for holding temporary values
 *  @ref_base:  the reference base
 *  @call:      filled with the annotations
 *
 *  Combines calls across the various samples being studied
 *  1. For each allele at each base across all samples the quality is summed so
 *     you end up with a set of quality sums for each allele present 2. The quality
 *     sums are sorted.
 *  3. Using the sorted quality sums we now create the allele ordering array
 *     A\subN. This is done by doing the following:
 *     a) If the reference allele is known it always comes first, otherwise N
 *        comes first.
 *     b) Then the rest of the alleles are output in descending order of quality
 *        sum (which we already know the qsum array was sorted).  Any allelles with
 *        qsum 0 will be excluded.
 *  4. Using the allele ordering array we create the genotype ordering array.
 *     In the worst case with an unknown reference this will be:  A0/A0 A1/A0 A1/A1
 *     A2/A0 A2/A1 A2/A2 A3/A0 A3/A1 A3/A2 A3/A3 A4/A0 A4/A1 A4/A2 A4/A3 A4/A4
 *  5. The genotype ordering array is then used to extract data from the error
 *     model 5*5 matrix and is used to produce a Phread likelihood array for each
 *     sample.
 */
int bcf_call_combine(int n, const bcf_callret1_t *calls, bcf_callaux_t *bca, int ref_base /*4-bit*/, bcf_call_t *call)
{
    int ref4, i, j;
    float qsum[5] = {0,0,0,0,0};
    if (ref_base >= 0) {
        call->ori_ref = ref4 = seq_nt16_int[ref_base];
        if (ref4 > 4) ref4 = 4;
    } else call->ori_ref = -1, ref4 = 0;

    // calculate qsum, this is done by summing normalized qsum across all samples,
    // to account for differences in coverage
    for (i = 0; i < n; ++i)
    {
        float sum = 0;
        for (j = 0; j < 4; ++j) sum += calls[i].qsum[j];
        if ( sum )
            for (j = 0; j < 4; j++) qsum[j] += calls[i].qsum[j] / sum;
    }

    // sort qsum in ascending order (insertion sort)
    float *ptr[5], *tmp;
    for (i=0; i<5; i++) ptr[i] = &qsum[i];
    for (i=1; i<4; i++)
        for (j=i; j>0 && *ptr[j] < *ptr[j-1]; j--)
            tmp = ptr[j], ptr[j] = ptr[j-1], ptr[j-1] = tmp;

    // Set the reference allele and alternative allele(s)
    for (i=0; i<5; i++) call->a[i] = -1;
    for (i=0; i<5; i++) call->qsum[i] = 0;
    call->unseen = -1;
    call->a[0] = ref4;
    for (i=3, j=1; i>=0; i--)   // i: alleles sorted by QS; j, a[j]: output allele ordering
    {
        int ipos = ptr[i] - qsum;   // position in sorted qsum array
        if ( ipos==ref4 )
            call->qsum[0] = qsum[ipos];    // REF's qsum
        else
        {
            if ( !qsum[ipos] ) break;       // qsum is 0, this and consequent alleles are not seen in the pileup
            call->qsum[j] = qsum[ipos];
            call->a[j++]  = ipos;
        }
    }
    if (ref_base >= 0)
    {
        // for SNPs, find the "unseen" base
        if (((ref4 < 4 && j < 4) || (ref4 == 4 && j < 5)) && i >= 0)
            call->unseen = j, call->a[j++] = ptr[i] - qsum;
        call->n_alleles = j;
    }
    else
    {
        call->n_alleles = j;
        if (call->n_alleles == 1) return -1; // no reliable supporting read. stop doing anything
    }
    /*
     * Set the phread likelihood array (call->PL) This array is 15 entries long
     * for each sample because that is size of an upper or lower triangle of a
     * worst case 5x5 matrix of possible genotypes. This worst case matrix will
     * occur when all 4 possible alleles are present and the reference allele
     * is unknown.  The sides of the matrix will correspond to the reference
     * allele (if known) followed by the alleles present in descending order of
     * quality sum
     */
    {
        int x, g[15], z;
        double sum_min = 0.;
        x = call->n_alleles * (call->n_alleles + 1) / 2;
        // get the possible genotypes
        // this is done by creating an ordered list of locations g for call (allele a, allele b) in the genotype likelihood matrix
        for (i = z = 0; i < call->n_alleles; ++i) {
            for (j = 0; j <= i; ++j) {
                g[z++] = call->a[j] * 5 + call->a[i];
            }
        }
        // for each sample calculate the PL
        for (i = 0; i < n; ++i)
        {
            int32_t *PL = call->PL + x * i;
            const bcf_callret1_t *r = calls + i;
            float min = FLT_MAX;
            for (j = 0; j < x; ++j) {
                if (min > r->p[g[j]]) min = r->p[g[j]];
            }
            sum_min += min;
            for (j = 0; j < x; ++j) {
                int y;
                y = (int)(r->p[g[j]] - min + .499);
                if (y > 255) y = 255;
                PL[j] = y;
            }
        }
        if ( call->DP4 )
        {
            for (i=0; i<n; i++)
            {
                call->DP4[4*i]   = calls[i].anno[0];
                call->DP4[4*i+1] = calls[i].anno[1];
                call->DP4[4*i+2] = calls[i].anno[2];
                call->DP4[4*i+3] = calls[i].anno[3];
            }
        }
        if ( call->ADF )
        {
            assert( call->n_alleles<=B2B_MAX_ALLELES );   // this is always true for SNPs and so far for indels as well

            // reorder ADR,ADF to match the allele ordering at this site
            int32_t tmp[B2B_MAX_ALLELES];
            int32_t *adr = call->ADR + B2B_MAX_ALLELES, *adr_out = call->ADR + B2B_MAX_ALLELES;
            int32_t *adf = call->ADF + B2B_MAX_ALLELES, *adf_out = call->ADF + B2B_MAX_ALLELES;
            int32_t *adr_tot = call->ADR;   // the first bin stores total counts per site
            int32_t *adf_tot = call->ADF;
            for (i=0; i<n; i++)
            {
                for (j=0; j<call->n_alleles; j++)
                {
                    tmp[j] = adr[ call->a[j] ];
                    adr_tot[j] += tmp[j];
                }
                for (j=0; j<call->n_alleles; j++) adr_out[j] = tmp[j];
                for (j=0; j<call->n_alleles; j++)
                {
                    tmp[j] = adf[ call->a[j] ];
                    adf_tot[j] += tmp[j];
                }
                for (j=0; j<call->n_alleles; j++) adf_out[j] = tmp[j];
                adf_out += call->n_alleles;
                adr_out += call->n_alleles;
                adr += B2B_MAX_ALLELES;
                adf += B2B_MAX_ALLELES;
            }
        }

//      if (ref_base < 0) fprintf(stderr, "%d,%d,%f,%d\n", call->n_alleles, x, sum_min, call->unseen);
        call->shift = (int)(sum_min + .499);
    }
    // combine annotations
    memset(call->anno, 0, 16 * sizeof(double));
    call->ori_depth = 0;
    call->depth     = 0;
    call->mq0       = 0;
    for (i = 0; i < n; ++i) {
        call->depth += calls[i].anno[0] + calls[i].anno[1] + calls[i].anno[2] + calls[i].anno[3];
        call->ori_depth += calls[i].ori_depth;
        call->mq0 += calls[i].mq0;
        for (j = 0; j < 16; ++j) call->anno[j] += calls[i].anno[j];
    }

    calc_SegBias(calls, call);

    // calc_chisq_bias("XPOS", call->bcf_hdr->id[BCF_DT_CTG][call->tid].key, call->pos, bca->ref_pos, bca->alt_pos, bca->npos);
    // calc_chisq_bias("XMQ", call->bcf_hdr->id[BCF_DT_CTG][call->tid].key, call->pos, bca->ref_mq, bca->alt_mq, bca->nqual);
    // calc_chisq_bias("XBQ", call->bcf_hdr->id[BCF_DT_CTG][call->tid].key, call->pos, bca->ref_bq, bca->alt_bq, bca->nqual);

    call->mwu_pos = calc_mwu_bias(bca->ref_pos, bca->alt_pos, bca->npos);
    call->mwu_mq  = calc_mwu_bias(bca->ref_mq,  bca->alt_mq,  bca->nqual);
    call->mwu_bq  = calc_mwu_bias(bca->ref_bq,  bca->alt_bq,  bca->nqual);
    call->mwu_mqs = calc_mwu_bias(bca->fwd_mqs, bca->rev_mqs, bca->nqual);

#if CDF_MWU_TESTS
    call->mwu_pos_cdf = calc_mwu_bias_cdf(bca->ref_pos, bca->alt_pos, bca->npos);
    call->mwu_mq_cdf  = calc_mwu_bias_cdf(bca->ref_mq,  bca->alt_mq,  bca->nqual);
    call->mwu_bq_cdf  = calc_mwu_bias_cdf(bca->ref_bq,  bca->alt_bq,  bca->nqual);
    call->mwu_mqs_cdf = calc_mwu_bias_cdf(bca->fwd_mqs, bca->rev_mqs, bca->nqual);
#endif

    call->vdb = calc_vdb(bca->alt_pos, bca->npos);

    return 0;
}

int bcf_call2bcf(bcf_call_t *bc, bcf1_t *rec, bcf_callret1_t *bcr, int fmt_flag, const bcf_callaux_t *bca, const char *ref)
{
    extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);
    int i, j, nals = 1;

    bcf_hdr_t *hdr = bc->bcf_hdr;
    rec->rid  = bc->tid;
    rec->pos  = bc->pos;
    rec->qual = 0;

    bc->tmp.l = 0;
    if (bc->ori_ref < 0)    // indel
    {
        // REF
        kputc(ref[bc->pos], &bc->tmp);
        for (j = 0; j < bca->indelreg; ++j) kputc(ref[bc->pos+1+j], &bc->tmp);

        // ALT
        for (i=1; i<4; i++)
        {
            if (bc->a[i] < 0) break;
            kputc(',', &bc->tmp); kputc(ref[bc->pos], &bc->tmp);

            if (bca->indel_types[bc->a[i]] < 0) { // deletion
                for (j = -bca->indel_types[bc->a[i]]; j < bca->indelreg; ++j)
                    kputc(ref[bc->pos+1+j], &bc->tmp);
            } else { // insertion; cannot be a reference unless a bug
                char *inscns = &bca->inscns[bc->a[i] * bca->maxins];
                for (j = 0; j < bca->indel_types[bc->a[i]]; ++j)
                    kputc("ACGTN"[(int)inscns[j]], &bc->tmp);
                for (j = 0; j < bca->indelreg; ++j) kputc(ref[bc->pos+1+j], &bc->tmp);
            }
            nals++;
        }
    }
    else    // SNP
    {
        kputc("ACGTN"[bc->ori_ref], &bc->tmp);
        for (i=1; i<5; i++)
        {
            if (bc->a[i] < 0) break;
            kputc(',', &bc->tmp);
            if ( bc->unseen==i ) kputs("<*>", &bc->tmp);
            else kputc("ACGT"[bc->a[i]], &bc->tmp);
            nals++;
        }
    }
    bcf_update_alleles_str(hdr, rec, bc->tmp.s);

    bc->tmp.l = 0;

    // INFO
    if (bc->ori_ref < 0)
    {
        bcf_update_info_flag(hdr, rec, "INDEL", NULL, 1);
        bcf_update_info_int32(hdr, rec, "IDV", &bca->max_support, 1);
        bcf_update_info_float(hdr, rec, "IMF", &bca->max_frac, 1);
    }
    bcf_update_info_int32(hdr, rec, "DP", &bc->ori_depth, 1);
    if ( fmt_flag&B2B_INFO_ADF )
        bcf_update_info_int32(hdr, rec, "ADF", bc->ADF, rec->n_allele);
    if ( fmt_flag&B2B_INFO_ADR )
        bcf_update_info_int32(hdr, rec, "ADR", bc->ADR, rec->n_allele);
    if ( fmt_flag&(B2B_INFO_AD|B2B_INFO_DPR) )
    {
        for (i=0; i<rec->n_allele; i++) bc->ADF[i] += bc->ADR[i];
        if ( fmt_flag&B2B_INFO_AD )
            bcf_update_info_int32(hdr, rec, "AD", bc->ADF, rec->n_allele);
        if ( fmt_flag&B2B_INFO_DPR )
            bcf_update_info_int32(hdr, rec, "DPR", bc->ADF, rec->n_allele);
    }

    float tmpf[16];
    for (i=0; i<16; i++) tmpf[i] = bc->anno[i];
    bcf_update_info_float(hdr, rec, "I16", tmpf, 16);
    bcf_update_info_float(hdr, rec, "QS", bc->qsum, nals);

    if ( bc->vdb != HUGE_VAL )      bcf_update_info_float(hdr, rec, "VDB", &bc->vdb, 1);
    if ( bc->seg_bias != HUGE_VAL ) bcf_update_info_float(hdr, rec, "SGB", &bc->seg_bias, 1);
    if ( bc->mwu_pos != HUGE_VAL )  bcf_update_info_float(hdr, rec, "RPB", &bc->mwu_pos, 1);
    if ( bc->mwu_mq != HUGE_VAL )   bcf_update_info_float(hdr, rec, "MQB", &bc->mwu_mq, 1);
    if ( bc->mwu_mqs != HUGE_VAL )  bcf_update_info_float(hdr, rec, "MQSB", &bc->mwu_mqs, 1);
    if ( bc->mwu_bq != HUGE_VAL )   bcf_update_info_float(hdr, rec, "BQB", &bc->mwu_bq, 1);
#if CDF_MWU_TESTS
    if ( bc->mwu_pos_cdf != HUGE_VAL )  bcf_update_info_float(hdr, rec, "RPB2", &bc->mwu_pos_cdf, 1);
    if ( bc->mwu_mq_cdf != HUGE_VAL )   bcf_update_info_float(hdr, rec, "MQB2", &bc->mwu_mq_cdf, 1);
    if ( bc->mwu_mqs_cdf != HUGE_VAL )  bcf_update_info_float(hdr, rec, "MQSB2", &bc->mwu_mqs_cdf, 1);
    if ( bc->mwu_bq_cdf != HUGE_VAL )   bcf_update_info_float(hdr, rec, "BQB2", &bc->mwu_bq_cdf, 1);
#endif
    tmpf[0] = bc->ori_depth ? (float)bc->mq0/bc->ori_depth : 0;
    bcf_update_info_float(hdr, rec, "MQ0F", tmpf, 1);

    // FORMAT
    rec->n_sample = bc->n;
    bcf_update_format_int32(hdr, rec, "PL", bc->PL, nals*(nals+1)/2 * rec->n_sample);
    if ( fmt_flag&B2B_FMT_DP )
    {
        int32_t *ptr = (int32_t*) bc->fmt_arr;
        for (i=0; i<bc->n; i++)
            ptr[i] = bc->DP4[4*i] + bc->DP4[4*i+1] + bc->DP4[4*i+2] + bc->DP4[4*i+3];
        bcf_update_format_int32(hdr, rec, "DP", bc->fmt_arr, rec->n_sample);
    }
    if ( fmt_flag&B2B_FMT_DV )
    {
        int32_t *ptr = (int32_t*) bc->fmt_arr;
        for (i=0; i<bc->n; i++)
            ptr[i] = bc->DP4[4*i+2] + bc->DP4[4*i+3];
        bcf_update_format_int32(hdr, rec, "DV", bc->fmt_arr, rec->n_sample);
    }
    if ( fmt_flag&B2B_FMT_SP )
    {
        int32_t *ptr = (int32_t*) bc->fmt_arr;
        for (i=0; i<bc->n; i++)
        {
            int fwd_ref = bc->DP4[4*i], rev_ref =  bc->DP4[4*i+1], fwd_alt = bc->DP4[4*i+2], rev_alt = bc->DP4[4*i+3];
            if ( fwd_ref+rev_ref<2 || fwd_alt+rev_alt<2 || fwd_ref+fwd_alt<2 || rev_ref+rev_alt<2 )
                ptr[i] = 0;
            else
            {
                double left, right, two;
                kt_fisher_exact(fwd_ref, rev_ref, fwd_alt, rev_alt, &left, &right, &two);
                int32_t x = (int)(-4.343 * log(two) + .499);
                if (x > 255) x = 255;
                ptr[i] = x;
            }
        }
        bcf_update_format_int32(hdr, rec, "SP", bc->fmt_arr, rec->n_sample);
    }
    if ( fmt_flag&B2B_FMT_DP4 )
        bcf_update_format_int32(hdr, rec, "DP4", bc->DP4, rec->n_sample*4);
    if ( fmt_flag&B2B_FMT_ADF )
        bcf_update_format_int32(hdr, rec, "ADF", bc->ADF+B2B_MAX_ALLELES, rec->n_sample*rec->n_allele);
    if ( fmt_flag&B2B_FMT_ADR )
        bcf_update_format_int32(hdr, rec, "ADR", bc->ADR+B2B_MAX_ALLELES, rec->n_sample*rec->n_allele);
    if ( fmt_flag&(B2B_FMT_AD|B2B_FMT_DPR) )
    {
        for (i=0; i<rec->n_sample*rec->n_allele; i++) bc->ADF[B2B_MAX_ALLELES+i] += bc->ADR[B2B_MAX_ALLELES+i];
        if ( fmt_flag&B2B_FMT_AD )
            bcf_update_format_int32(hdr, rec, "AD", bc->ADF+B2B_MAX_ALLELES, rec->n_sample*rec->n_allele);
        if ( fmt_flag&B2B_FMT_DPR )
            bcf_update_format_int32(hdr, rec, "DPR", bc->ADF+B2B_MAX_ALLELES, rec->n_sample*rec->n_allele);
    }

    return 0;
}
