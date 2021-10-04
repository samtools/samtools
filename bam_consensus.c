// FIXME: also use strand to spot possible basecalling errors.
//        Specifically het calls where mods are predominantly on one
//        strand.  So maybe require + and - calls and check concordance
//        before calling a het as confident.  (Still call, but low qual?)

// FIXME: add bedmethyl support


// Eg 50T 20A seems T/A het,
// but 30T+ 20T- 18A+ 2A- seems like a consistent A miscall on one strand
// only, while T is spread evenly across both strands.

/*  bam_consensus.c -- consensus subcommand.

    Copyright (C) 1998-2001,2003 Medical Research Council (Gap4/5 source)
    Copyright (C) 2003-2005,2007-2021 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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

/*
 * The origins of this algorithm come from the Staden Package, initially
 * the Gap4 consensus algorithm, revised for Gap5, and substantially
 * rewritten again for the qual-lossy Crumble compression tool which
 * is where this algorithm has forked from.
 */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>

#include <htslib/sam.h>

#include "samtools.h"
#include "sam_opts.h"

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

// Minimum cutoff for storing mod data; => at least 10% chance
#define MOD_CUTOFF 0.46

enum format {
    FASTQ,
    FASTA,
    PILEUP
};

typedef unsigned char uc;

typedef struct {
    samFile *fp;
    sam_hdr_t *h;
    hts_idx_t *idx;
    hts_itr_t *iter;
    char *reg;
    int use_qual;
    int use_mqual;
    int min_depth;
    double call_fract;
    double het_fract;
    int gap5;
    enum format fmt;
    int cons_cutoff;
    int ambig;
    int line_len;
    int mods;
    int default_qual;
    int mod_cutoff;  // phred-style cutoff for calling modified base
    double mod_prob; // used in computing mod quality; Avg P(call) per event.
    double mod_yes;  // >=Q threshold for asserting mod call to be yes
    double mod_no;   // <=Q threshold for asserting mod call to be no
    int het_only;
    int all_bases;
    int show_del;
    int show_ins;
} consensus_opts;

static int readaln(void *data, bam1_t *b) {
    consensus_opts *dat = (consensus_opts *)data;
    if (dat->iter)
        return sam_itr_next(dat->fp, dat->iter, b);
    else
        return sam_read1(dat->fp, dat->h, b);
}

// Initialise and destroy the base modifier state data. This is called
// as each new read is added or removed from the pileups.
int pileup_cd_create(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    hts_base_mod_state *m = hts_base_mod_state_alloc();
    bam_parse_basemod(b, m);
    cd->p = m;
    return 0;
}

int pileup_cd_destroy(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    hts_base_mod_state_free(cd->p);
    return 0;
}

/* --------------------------------------------------------------------------
 * A bayesian consensus algorithm that analyses the data to work out
 * which hypothesis of pure A/C/G/T/absent and all combinations of two
 * such bases meets the observations.
 *
 * This has its origins in Gap4 (homozygous) -> Gap5 (heterozygous)
 * -> Crumble (tidied up to use htslib's pileup) -> here.
 *
 */

#define CONS_DISCREP    4
#define CONS_ALL        15

#define CONS_MQUAL      16

// Make these a parameter, but capping mqual seems wise as it's
// not so robust especially at the low end when real SNPs exist,
// but reduce mqual for many algorithms.
#define MIN_MQUAL       10
#define MAX_MQUAL       50

typedef struct {
    /* the most likely base call - we never call N here */
    /* A=0, C=1, G=2, T=3, *=4 */
    int call;

    /* The most likely heterozygous base call */
    /* Use "ACGT*"[het / 5] vs "ACGT*"[het % 5] for the combination */
    int het_call;

    /* Log-odds for het_call */
    int het_logodd;

    /* Single phred style call */
    int phred;

    /* Sequence depth */
    int depth;

    /* Discrepancy search score */
    float discrep;

    /* Base modifications, for top (0) and bottom (1) strands */
    int mod_call[2];  // +ve for char, -ve for ChEBI ID
    int mod_qual[2];  // phred scale
} consensus_t;

static char *seq_nt16_str_pad = "=ACMGRSVTWYHKDBN=acmgrsvtwyhkdbn";

/*
 * Produce the insertion portion of the consensus sequence.
 * Returns 0 on success, filling out ins_ks and ins_ks,
 *        -1 on failure
 */
static int insertion_consensus(const bam_pileup1_t *plp, int nplp,
                               consensus_opts *opts, int nins, int ins_len,
                               kstring_t *ins_ks, kstring_t *inq_ks) {
    // Map "seqi" nt16 to A,C,G,T compatibility with weights on pure bases.
    // where seqi is A | (C<<1) | (G<<2) | (T<<3)
    //                        * A C M  G R S V  T W Y H  K D B N
    static int seqi2A[16] = { 0,8,0,4, 0,4,0,2, 0,4,0,2, 0,2,0,1 };
    static int seqi2C[16] = { 0,0,8,4, 0,0,4,2, 0,0,4,2, 0,0,2,1 };
    static int seqi2G[16] = { 0,0,0,0, 8,4,4,1, 0,0,0,0, 4,2,2,1 };
    static int seqi2T[16] = { 0,0,0,0, 0,0,0,0, 8,4,4,2, 8,2,2,1 };

    int ins_[128][16], (*ins)[16] = ins_;
    int i;

    if (ins_len > 128)
        if (!(ins = calloc(ins_len * 16, sizeof(int))))
            return -1;

    // Accumulate into ins[][] arrays
    kstring_t is = {0,0};
    for (i = 0; i < nplp; i++) {
        const bam_pileup1_t *p = plp+i;

        // Slow method for now is bam_plp_insertion and parse back.
        // NB no support for opts->use_qual yet here.  We either need
        // a better bam_plp_insertion or inline it here with added
        // quality value processing.
        int il = bam_plp_insertion(p, &is, NULL), j;
        if (il != ins_len)
            continue;

        for (j = 0; j < il; j++) {
            if (is.s[j] == '*') {
                ins[j][0] += 8;
            } else {
                int Q = seqi2A[seq_nt16_table[(uc)is.s[j]]];
                ins[j][1] += Q;
                Q = seqi2C[seq_nt16_table[(uc)is.s[j]]];
                ins[j][2] += Q;
                Q = seqi2G[seq_nt16_table[(uc)is.s[j]]];
                ins[j][4] += Q;
                Q = seqi2T[seq_nt16_table[(uc)is.s[j]]];
                ins[j][8] += Q;
            }
        }
    }

    // Call insertion consensus from ins[][]
    int j;
    ins_ks->l = 0;
    for (j = 0; j < ins_len; j++) {
        int k, b1 = 0, b2 = 0, b3 = 0, m1 = 0, m2 = 0, m3 = 0;
        for (k = 0; k < 16; k++) {
            if (m1 < ins[j][k])
                m3 = m2, m2 = m1, m1 = ins[j][k], b2 = b1, b1 = k;
            else if (m2 < ins[j][k])
                m3 = m2, m2 = ins[j][k], b2 = k;
            else if (m3 < ins[j][k])
                m3 = ins[j][k], b3 = k;
        }
        if (b1+b2+b3 != 0) {
            int call = b1 ? b1 : 16;
            int score = m1;
            if (b2 != 0 && opts->ambig && m2 >= opts->het_fract * m1) {
                call |= b2 ? b2 : 16;
                score += m2;
            }
            if (b3 != 0 && opts->ambig && m3 >= opts->het_fract * m1) {
                call |= b3 ? b3 : 16;
                score += m3;
            }

            char b = opts->call_fract-DBL_EPSILON >= score / (nins*8.0)
                ? 'N' : seq_nt16_str_pad[call];
            kputc(b, ins_ks);
            int qual = 100 * m1 / nplp;
            qual = qual + '!' > '~' ? '~' : qual + '!';
            kputc(qual, inq_ks);
        }
    }

    if (ins != ins_)
        free(ins);

    return 0;
}

#define P_HET 1e-6

#define LOG10            2.30258509299404568401
#define TENOVERLOG10     4.34294481903251827652
#define TENLOG2OVERLOG10 3.0103


// FIXME: not known at present.
// Best to externalise these with parameters?  Ideally would be
// param per read-group though.

/* Sequencing technologies for seq_t.seq_tech; 5 bits, so max=31 */
#define STECH_UNKNOWN    0
#define STECH_SANGER     1
#define STECH_SOLEXA     2
#define STECH_SOLID      3
#define STECH_454        4
#define STECH_HELICOS    5
#define STECH_IONTORRENT 6
#define STECH_PACBIO     7
#define STECH_ONT        8
#define STECH_LAST       8 // highest value

// Undercall rates govern the alignment pad (deletion) vs missing observation
// (undercall).  Single molecule techniques may miss data, but we don't want
// to start labelling lots of places as being base/gap het calls.
//
// NB: Figures are pure guesses, and also unused! (forced to illumina atm)
double tech_undercall[] = {
    1.00, // unknown
    1.00, // sanger
    1.00, // solexa/illumina
    1.00, // solid
    1.00, // 454
    1.00, // helicos
    1.00, // iontorrent
    1.00, // pacbio
    1.20, // ont
};

#ifdef __GNUC__
#define ALIGNED(x) __attribute((aligned(x)))
#else
#define ALIGNED(x)
#endif

static double prior[25]    ALIGNED(16);  /* Sum to 1.0 */
static double lprior15[15] ALIGNED(16);  /* 15 combinations of {ACGT*} */

/* Precomputed matrices for the consensus algorithm */
static double pMM[9][101] ALIGNED(16);
static double p__[9][101] ALIGNED(16);
static double p_M[9][101] ALIGNED(16);
static double po_[9][101] ALIGNED(16);
static double poM[9][101] ALIGNED(16);
static double poo[9][101] ALIGNED(16);
static double puu[9][101] ALIGNED(16);
static double pum[9][101] ALIGNED(16);
static double pmm[9][101] ALIGNED(16);

static double e_tab_a[1002]  ALIGNED(16);
static double *e_tab = &e_tab_a[500];
static double e_tab2_a[1002] ALIGNED(16);
static double *e_tab2 = &e_tab2_a[500];
static double e_log[501]     ALIGNED(16);

/*
 * Lots of confusing matrix terms here, so some definitions will help.
 *
 * M = match base
 * m = match pad
 * _ = mismatch
 * o = overcall
 * u = undercall
 *
 * We need to distinguish between homozygous columns and heterozygous columns,
 * done using a flat prior.  This is implemented by treating every observation
 * as coming from one of two alleles, giving us a 2D matrix of possibilities
 * (the hypotheses) for each and every call (the observation).
 *
 * So pMM[] is the chance that given a call 'x' that it came from the
 * x/x allele combination.  Similarly p_o[] is the chance that call
 * 'x' came from a mismatch (non-x) / overcall (consensus=*) combination.
 *
 * Examples with observation (call) C and * follows
 *
 *  C | A  C  G  T  *          * | A  C  G  T  *
 *  -----------------          -----------------
 *  A | __ _M __ __ o_         A | uu uu uu uu um
 *  C | _M MM _M _M oM         C | uu uu uu uu um
 *  G | __ _M __ __ o_         G | uu uu uu uu um
 *  T | __ _M __ __ o_         T | uu uu uu uu um
 *  * | o_ oM o_ o_ oo         * | um um um um mm
 *
 * In calculation terms, the _M is half __ and half MM, similarly o_ and um.
 *
 * Relative weights of substitution vs overcall vs undercall are governed on a
 * per base basis using the P_OVER and P_UNDER scores (subst is
 * 1-P_OVER-P_UNDER).
 *
 * The heterozygosity weight though is a per column calculation as we're
 * trying to model whether the column is pure or mixed. Hence this is done
 * once via a prior and has no affect on the individual matrix cells.
 */

static void consensus_init(double p_het) {
    int i, t;

    for (i = -500; i <= 500; i++)
        e_tab[i] = exp(i);
    for (i = -500; i <= 500; i++)
        e_tab2[i] = exp(i/10.);
    for (i = 0; i <= 500; i++)
        e_log[i] = log(i);

    // Heterozygous locations
    for (i = 0; i < 25; i++)
        prior[i] = p_het / 20;
    prior[0] = prior[6] = prior[12] = prior[18] = prior[24] = (1-p_het)/5;

    lprior15[0]  = log(prior[0]);
    lprior15[1]  = log(prior[1]*2);
    lprior15[2]  = log(prior[2]*2);
    lprior15[3]  = log(prior[3]*2);
    lprior15[4]  = log(prior[4]*2);
    lprior15[5]  = log(prior[6]);
    lprior15[6]  = log(prior[7]*2);
    lprior15[7]  = log(prior[8]*2);
    lprior15[8]  = log(prior[9]*2);
    lprior15[9]  = log(prior[12]);
    lprior15[10] = log(prior[13]*2);
    lprior15[11] = log(prior[14]*2);
    lprior15[12] = log(prior[18]);
    lprior15[13] = log(prior[19]*2);
    lprior15[14] = log(prior[24]);


    // Rewrite as new form
    for (t = STECH_UNKNOWN; t <= STECH_LAST; t++) {
        for (i = 1; i < 101; i++) {
            double prob = 1 - pow(10, -i / 10.0);

            // May want to multiply all these by 5 so pMM[i] becomes close
            // to -0 for most data. This makes the sums increment very slowly,
            // keeping bit precision in the accumulator.
            pMM[t][i] = log(prob/5);
            p__[t][i] = log((1-prob)/20);
            p_M[t][i] = log((exp(pMM[t][i]) + exp(p__[t][i]))/2);

            puu[t][i] = p__[t][i];

            poM[t][i] = p_M[t][i] *= tech_undercall[t];
            po_[t][i] = p__[t][i] *= tech_undercall[t];
            poo[t][i] = p__[t][i] *= tech_undercall[t];
            pum[t][i] = p_M[t][i] *= tech_undercall[t];
            pmm[t][i] = pMM[t][i] *= tech_undercall[t];
        }

        pMM[t][0] = pMM[t][1];
        p__[t][0] = p__[t][1];
        p_M[t][0] = p_M[t][1];

        pmm[t][0] = pmm[t][1];
        poo[t][0] = poo[t][1];
        po_[t][0] = po_[t][1];
        poM[t][0] = poM[t][1];
        puu[t][0] = puu[t][1];
        pum[t][0] = pum[t][1];
    }
}

static inline double fast_exp(double y) {
    if (y >= -50 && y <= 50)
        return e_tab2[(int)(y*10)];

    if (y < -500)
        y = -500;
    if (y > 500)
        y = 500;

    return e_tab[(int)y];
}

/* Taylor (deg 3) implementation of the log:
 * http://www.flipcode.com/cgi-bin/fcarticles.cgi?show=63828
 */
static inline double fast_log2(double val)
{
   register int64_t *const     exp_ptr = ((int64_t*)&val);
   register int64_t            x = *exp_ptr;
   register const int      log_2 = ((x >> 52) & 2047) - 1024;
   x &= ~(2047LL << 52);
   x += 1023LL << 52;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;

   return val + log_2;
}

#define ph_log(x) (-TENLOG2OVERLOG10*fast_log2((x)))


/*
 * In mammals, nearly all C mods are in CG sequences (CpG), with Cs on
 * both strands being modified: https://en.wikipedia.org/wiki/DNA_methylation
 *
 * Initial nanoporetech/megalodon appears to follow this, only ever
 * emitting CG to mG calls.  Furthermore it looks like 100% of such
 * occurences are called, but obviously the probability values are
 * mostly low indicating it believes they're unlikely to be true.
 *
 * This poses a problem though, as half the time CG is CG(fwd) and
 * GC(rev) meaning only half of them can be used for determining
 * probability of the mod occurring.
 *
 * The solution is to filter mods by strand. So if it's C[+m<qual>]
 * then we count the mods against other fwd strand C only.
 *
 * The cases where it may miss a CG and not emit a call are when the
 * read starts or ends between it.  Ie C$ or g^? in mpileup terms.
 */

/*
 * As per calculate_consensus_bit_het but for a single pileup column.
 */
static
int calculate_consensus_gap5(int pos, int flags,
                             const bam_pileup1_t *p,
                             int np,
                             consensus_opts *opts,
                             consensus_t *cons,
                             int default_qual,
                             kstring_t *mod_ks,
                             double mod_prob,
                             double mod_yes_threshold,
                             double mod_no_threshold,
                             kstring_t *ins_ks,
                             kstring_t *inq_ks) {
    int i, j;
    static int init_done =0;
    static double q2p[101], mqual_pow[256];
    double min_e_exp = DBL_MIN_EXP * log(2) + 1;

    double S[15] ALIGNED(16) = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double sumsC[6] = {0,0,0,0,0,0};
    int depth = 0;

    /* Base mods basic strategy:
     *
     * Two things to consider; modified or not?
     * If modified, then what is the most likely modification?
     * We also do this per strand.
     */
    int mods_present[2][257];  // list ending in 0
    int mod_canonical[2][256]; // which canonical base type to use
    double mod[2][256][3];     // [strand][type][called y/n/unknown]
    int unmod[2][15];          // [strand][base]

    if (mod_ks) {
        mods_present[0][0] = 0;
        mods_present[1][0] = 0;
        mod[0][0][0] = 0;
        mod[0][0][1] = 0;
        mod[0][0][2] = 0;
        mod[1][0][0] = 0;
        mod[1][0][1] = 0;
        mod[1][0][2] = 0;

        memset(unmod, 0, 2*15*sizeof(int));
    }

    /* Map the 15 possible combinations to 1-base or 2-base encodings */
    static int map_sing[15] ALIGNED(16) =
        {0, 5, 5, 5, 5,
            1, 5, 5, 5,
               2, 5, 5,
                  3, 5,
                     4};
    static int map_het[15] ALIGNED(16) =
        {0,  1,  2,  3,  4,
             6,  7,  8,  9,
                12, 13, 14,
                    18, 19,
                        24};

    if (!init_done) {
        init_done = 1;
        consensus_init(P_HET);

        for (i = 0; i <= 100; i++) {
            q2p[i] = pow(10, -i/10.0);
        }

        for (i = 0; i < 255; i++) {
            //mqual_pow[i] = 1-pow(10, -(i+.01)/10.0);
            mqual_pow[i] = 1-pow(10, -(i*.9)/10.0);
            //mqual_pow[i] = 1-pow(10, -(i/3+.1)/10.0);
            //mqual_pow[i] = 1-pow(10, -(i/2+.05)/10.0);
        }
        // unknown mqual
        mqual_pow[255] = mqual_pow[10];
    }

    /* Initialise */
    int counts[6] = {0};

    /* Accumulate */
    // FIXME: also seed with unknown alleles so low coverage data is
    // less confident.
    int ins_len[10] = {0}, ins_len_freq[10] = {0}, nins = 0;
    int n;
    for (n = 0; n < np; n++) {
        if (p[n].is_refskip)
            continue;

        bam1_t *b = p[n].b;
        uint8_t base = bam_seqi(bam_get_seq(b), p[n].qpos);
        uint8_t *qual_arr = bam_get_qual(b);
        uint8_t qual = qual_arr[p[n].qpos];
        if (qual == 255 || (qual == 0 && *qual_arr == 255))
            qual = default_qual;
        const int stech = STECH_SOLEXA;

        // =ACM GRSV TWYH KDBN
        static int L[16] = {
            5,0,1,5, 2,5,5,5, 3,5,5,5, 5,5,5,5
        };

        // convert from sam base to acgt*n order.
        base = L[base];
        if (p[n].is_del) base = 4;

        double MM, __, _M, qe;

        // Correction for mapping quality.  Maybe speed up via lookups?
        // Cannot nullify mapping quality completely.  Lots of (true)
        // SNPs means low mapping quality.  (Ideally need to know
        // hamming distance to next best location.)

        if (p[n].indel > 0) {
            // FIXME: add flag&CONS_MQUAL support too
            nins++;
            int j;

            // Build insertion length distribution:
            // We fill out ins_len[] (lengths) and ins_len_freq[] and
            // maintain as sorted top 10, with ele 0 as most frequent.

            // Length may need adjusting due to pads
            int indel_len = p[n].indel;
            uint32_t *cig = bam_get_cigar(p[n].b);
            for (j = p[n].cigar_ind+1; j < p[n].b->core.n_cigar; j++) {
                switch (cig[j] & BAM_CIGAR_MASK) {
                case BAM_CPAD:
                    indel_len += cig[j] >> BAM_CIGAR_SHIFT;
                    break;
                case BAM_CINS:
                    break;
                default:
                    j = p[n].b->core.n_cigar; // terminate for loop
                }
            }

            for (j = 0; j < 10 && ins_len[j]; j++)
                if (ins_len[j] == indel_len)
                    break;
            if (j < 10) {
                ins_len_freq[j]++;
                ins_len[j] = indel_len;
            }
            while (j > 0 && ins_len_freq[j] > ins_len_freq[j-1]) {
                int tmp = ins_len[j-1];
                ins_len[j-1] = ins_len[j];
                ins_len[j] = tmp;
                tmp = ins_len_freq[j-1];
                ins_len_freq[j-1] = ins_len_freq[j];
                ins_len_freq[j] = tmp;
            }
        }

        if (flags & CONS_MQUAL) {
            double _p = 1-q2p[qual];
            int mqual = b->core.qual;
            if (mqual < MIN_MQUAL)
                mqual = MIN_MQUAL;
            if (mqual > MAX_MQUAL)
                mqual = MAX_MQUAL;
            //double _m = mqual_pow[mqual];
            double _m = 1-q2p[mqual];

            qual = ph_log(1-(_m * _p + (1 - _m)/4));
            //qual = ph_log(1-_p*_m);
        }

        /* Quality 0 should never be permitted as it breaks the math */
        if (qual < 1)
            qual = 1;

        __ = p__[stech][qual];
        MM = pMM[stech][qual] - __;
        _M = p_M[stech][qual] - __;

        if (flags & CONS_DISCREP) {
            qe = q2p[qual];
            sumsC[base] += 1 - qe;
        }

        counts[base]++;

        switch (base) {
        case 0:
            S[0] += MM; S[1 ]+= _M; S[2 ]+= _M; S[3 ]+= _M; S[4 ]+= _M;
            break;

        case 1:
            S[1 ]+= _M; S[5 ]+= MM; S[6 ]+= _M; S[7 ]+= _M; S[8 ]+= _M;
            break;

        case 2:
            S[2 ]+= _M; S[6 ]+= _M; S[9 ]+= MM; S[10]+= _M; S[11]+= _M;
            break;

        case 3:
            S[3 ]+= _M; S[7 ]+= _M; S[10]+= _M; S[12]+= MM; S[13]+= _M;
            break;

        case 4:
            S[4 ]+= _M; S[8 ]+= _M; S[11]+= _M; S[13]+= _M; S[14]+= MM;
            break;

        case 5: /* N => equal weight to all A,C,G,T but not a pad */
            S[0] += MM; S[1 ]+= MM; S[2 ]+= MM; S[3 ]+= MM; S[4 ]+= _M;
                        S[5 ]+= MM; S[6 ]+= MM; S[7 ]+= MM; S[8 ]+= _M;
                                    S[9 ]+= MM; S[10]+= MM; S[11]+= _M;
                                                S[12]+= MM; S[13]+= _M;
            break;
        }

        depth++;

        // We do clipping on our probabitilies so we can judge the number
        // of highly confidently modified bases against the number of
        // highly confidently unmodified bases.
        // Very high confidence calls are used as positive evidence. [0]
        // Very low confidence calls are used as negative evidence.  [1]
        // Inbetween ones get put into the "unknown" category.       [2]
        //
        // These levels can be adjusted so users get the option of
        // controlling whether negative calls count against the consensus or
        // are simply ignored (as part of a heterozygous site).
        if (mod_ks && !p[n].is_del && p[n].indel <= 0) {
            int strand = b->core.flag & BAM_FREVERSE ? 1 : 0;
            hts_base_mod_state *m = p[n].cd.p;
            hts_base_mod mods[256];
            int nm = bam_mods_at_qpos(b, p[n].qpos, m, mods, 256);

            if (nm) {
                int j, k;
                double p_tot = 0;
                for (j = 0; j < nm && j < 256; j++) {
                    double prob = mods[j].qual / 256.0;
                    prob += p_tot;
                    for (k = 0; k < 256 && mod[strand][k][0]; k++)
                        if (mods_present[strand][k] == mods[j].modified_base)
                            break;
                    if (!mods_present[strand][k]) {
                        mods_present[strand][k+1] = 0;
                        mods_present[strand][k]   = mods[j].modified_base;
                        mod_canonical[strand][k]  = mods[j].canonical_base;
                        mod[strand][k][0] = 0;
                        mod[strand][k][1] = 0;
                        mod[strand][k][2] = 0;
                    }

                    // SUM(mod[strand][k][0..2]) = depth of called.
                    // Plus unmod[strand][base] for depth of uncalled
                    // to get total freq.
                    if (j == 0) {
                        // Counts with most likely call only.
                        if (prob >= mod_yes_threshold)
                            mod[strand][k][0]++;
                        else if (prob <= mod_no_threshold)
                            mod[strand][k][1]++;
                        else
                            mod[strand][k][2]++;
                    }
                }
            } else {
                // Else do something with MISSING_MODP.
                // We don't know which mod types will be present, so
                // we just record these under their canonical base type
                // for now and add it in later.
                unmod[strand][base]++;
            }
        }
    }

    // Insertion seqeuence
    if (ins_len_freq[0] && nins >= np/2 && nins >= opts->min_depth) {
        if (insertion_consensus(p, np, opts, nins, ins_len[0],
                                ins_ks, inq_ks) < 0)
            return -1;
    }

    /* We've accumulated stats, so now we speculate on the consensus call */
    {
        double shift, max, max_het, norm[15];
        int call = 0, het_call = 0, ph;
        double tot1, tot2;

        /*
         * Scale numbers so the maximum score is 0. This shift is essentially
         * a multiplication in non-log scale to both numerator and denominator,
         * so it cancels out. We do this to avoid calling exp(-large_num) and
         * ending up with norm == 0 and hence a 0/0 error.
         *
         * Can also generate the base-call here too.
         */
        shift = -DBL_MAX;
        max = -DBL_MAX;
        max_het = -DBL_MAX;

        for (j = 0; j < 15; j++) {
            S[j] += lprior15[j];
            if (shift < S[j])
                shift = S[j];

            /* Only call pure AA, CC, GG, TT, ** for now */
            if (j != 0 && j != 5 && j != 9 && j != 12 && j != 14) {
                if (max_het < S[j]) {
                    max_het = S[j];
                    het_call = j;
                }
                continue;
            }

            if (max < S[j]) {
                max = S[j];
                call = j;
            }
        }

        /*
         * Shift and normalise.
         * If call is, say, b we want p = b/(a+b+c+...+n), but then we do
         * p/(1-p) later on and this has exceptions when p is very close
         * to 1.
         *
         * Hence we compute b/(a+b+c+...+n - b) and
         * rearrange (p/norm) / (1 - (p/norm)) to be p/norm2.
         */
        for (j = 0; j < 15; j++) {
            S[j] -= shift;
            double e = fast_exp(S[j]);
            S[j] = (S[j] > min_e_exp) ? e : DBL_MIN;
            norm[j] = 0;
        }

        tot1 = tot2 = 0;
        for (j = 0; j < 15; j++) {
            norm[j]    += tot1;
            norm[14-j] += tot2;
            tot1 += S[j];
            tot2 += S[14-j];
        }

        /* And store result */
        if (depth && depth != counts[5] /* all N */) {
            double m;

            cons->depth = depth;

            cons->call     = map_sing[call];
            if (norm[call] == 0) norm[call] = DBL_MIN;
            // approximation of phred for when S[call] ~= 1 and norm[call]
            // is small.  Otherwise we need the full calculation.
            if (S[call] == 1 && norm[call] < .01)
                ph = ph_log(norm[call]) + .5;
            else
                ph = ph_log(1-S[call]/(norm[call]+S[call])) + .5;
            cons->phred = ph < 0 ? 0 : ph;
            //cons->call_prob1 = norm[call]; // p = 1 - call_prob1

            cons->het_call = map_het[het_call];
            if (norm[het_call] == 0) norm[het_call] = DBL_MIN;
            ph = TENLOG2OVERLOG10 * (fast_log2(S[het_call]) - fast_log2(norm[het_call])) + .5;

            //printf("Call=%c%c %2d, hetcall = %c%c %2d\n",
            //       "ACGT*"[map_het[call]%5], "ACGTN"[map_het[call]/5],
            //       cons->phred,
            //       "ACGT*"[het_call%5], "ACGTN"[het_call/5],
            //       ph);
            cons->het_logodd = ph;
            //cons->het_prob_n = S[het_call]; // p = prob_n / prob_d
            //cons->het_prob_d = norm[het_call];

            /* Compute discrepancy score */
            if (flags & CONS_DISCREP) {
                m = sumsC[0]+sumsC[1]+sumsC[2]+sumsC[3]+sumsC[4];
                double c;
                if (cons->het_logodd > 0)
                    c = sumsC[cons->het_call%5] + sumsC[cons->het_call/5];
                else
                    c = sumsC[cons->call];
                cons->discrep = (m-c)/sqrt(m);
            }
        } else {
            cons->call = 4; /* N */
            cons->het_call = 0;
            cons->het_logodd = 0;
            cons->phred = 0;
            cons->depth = 0;
            cons->discrep = 0;
        }
    }

    cons->mod_qual[0] = 0;
    cons->mod_qual[1] = 0;

    /* We have the consensus call, so now check for base modifications. */
    if (mod_ks) {
        int s, k;
        for (s = 0; s < 2; s++) { // strand
            char *rc = s==0 ? "ACGT*" : "TGCA*";
            double best_qual = 0;
            for (k = 0; mods_present[s][k]; k++) {
                if (mod_canonical[s][k] != rc[cons->call] &&
                    mod_canonical[s][k] != 'N') {
                    // If the base modification is not compatible with the
                    // consensus call, then we deem the modification to also
                    // be incorrect.  (It'll have arrived due to
                    // multiple base types aligned to this location.)
                    continue;
                    //fprintf(stderr, "[consensus]: detected non-matching "
                    //        "mod %c->%c vs consensus %c at pos %d\n",
                    //        mod_canonical[s][k], mods_present[s][k],
                    //        rc[cons->call], pos+1);
                }
                // Treat unmod bases as bases with a very low call probability.
                mod[s][k][1] += unmod[s][cons->call];

                // N = freq(mod_yes) - freq(mod_no)
                // X = base prob of frequency occurance being correct
                // P = X^N / (X^N + (1-X)^N)
                double N = mod[s][k][0] - mod[s][k][1];
                double X = mod_prob;
                double pqual = pow(X, N) / (pow(X, N) + pow(1-X, N));

                pqual = -10 * log(1-pqual+DBL_MIN) / log(10);
                if (pqual < 0) pqual = 0;

                // FIXME: Apply some correction with mod[s][k][2].  Too
                // many with indeterministic confidence is probably an
                // indication of failure still.
                //
                // Eg I'd trust 10Y, 2N much more than 50Y, 0N and 500 unk.

                if (pqual >= MOD_CUTOFF && pqual > best_qual) {
                    best_qual = pqual;
                    cons->mod_call[s] = mods_present[s][k];
                    cons->mod_qual[s] = pqual+.5;
                }
            }
        }
    }

    return 0;
}


/* --------------------------------------------------------------------------
 * A simple summing algorithm, either pure base frequency, or by
 * weighting them according to their quality values.
 *
 * This is crude, but easy to understand and fits with several
 * standard pileup criteria (eg COG-UK / CLIMB Covid-19 seq project).
 *
 *
 * call1 / score1 / depth1 is the highest scoring allele.
 * call2 / score2 / depth2 is the second highest scoring allele.
 *
 * Het_fract:  score2/score1
 * Call_fract: score1 or score1+score2 over total score
 * Min_depth:  minimum total depth of utilised bases (depth1+depth2)
 * Min_score:  minimum total score of utilised bases (score1+score2)
 *
 * Eg het_fract 0.66, call_fract 0.75 and min_depth 10.
 * 11A, 2C, 2G (14 total depth) is A.
 * 9A, 2C, 2G  (12 total depth) is N as depth(A) < 10.
 * 11A, 5C, 5G (21 total depth) is N as 11/21 < 0.75 (call_fract)
 *
 *
 * 6A, 5G, 1C  (12 total depth) is AG het as depth(A)+depth(G) >= 10
 *                              and 5/6 >= 0.66 and 11/12 >= 0.75.
 *
 * 6A, 5G, 4C  (15 total depth) is N as (6+5)/15 < 0.75 (call_fract).
 *
 *
 * Note for the purpose of deletions, a base/del has an ambiguity
 * code of lower-case base (otherwise it is uppercase).
 */

static int calculate_consensus_simple(const bam_pileup1_t *plp, int nplp,
                                      consensus_opts *opts, int *qual,
                                      kstring_t *ins_ks, kstring_t *inq_ks) {
    int i, min_qual = 0;

    // Map "seqi" nt16 to A,C,G,T compatibility with weights on pure bases.
    // where seqi is A | (C<<1) | (G<<2) | (T<<3)
    //                        * A C M  G R S V  T W Y H  K D B N
    static int seqi2A[16] = { 0,8,0,4, 0,4,0,2, 0,4,0,2, 0,2,0,1 };
    static int seqi2C[16] = { 0,0,8,4, 0,0,4,2, 0,0,4,2, 0,0,2,1 };
    static int seqi2G[16] = { 0,0,0,0, 8,4,4,1, 0,0,0,0, 4,2,2,1 };
    static int seqi2T[16] = { 0,0,0,0, 0,0,0,0, 8,4,4,2, 8,2,2,1 };

    // Ignore ambiguous bases in seq for now, so we don't treat R, Y, etc
    // as part of one base and part another.  Based on BAM seqi values.
    // We also use freq[16] as "*" for gap.
    int freq[17] = {0};  // base frequency, aka depth
    int score[17] = {0}; // summation of base qualities
    int ins_len[10] = {0}, ins_len_freq[10] = {0}, nins = 0;

    // Accumulate
    for (i = 0; i < nplp; i++) {
        const bam_pileup1_t *p = plp+i;
        int q = bam_get_qual(p->b)[p->qpos];
        if (q < min_qual)
            // Should we still record these in freq[] somewhere so
            // we can use them in the fracts?
            // Difference between >= X% of high-qual bases calling Y
            // and >= X% of all bases are high-quality Y calls.
            continue;

        int b = p->is_del ? 16 : bam_seqi(bam_get_seq(p->b), p->qpos);

        if (p->indel > 0) {
            nins++;
            int j;

            // Build insertion length distribution:
            // We fill out ins_len[] (lengths) and ins_len_freq[] and
            // maintain as sorted top 10, with ele 0 as most frequent.

            // Length may need adjusting due to pads
            int indel_len = p->indel;
            uint32_t *cig = bam_get_cigar(p->b);
            for (j = p->cigar_ind+1; j < p->b->core.n_cigar; j++) {
                switch (cig[j] & BAM_CIGAR_MASK) {
                case BAM_CPAD:
                    indel_len += cig[j] >> BAM_CIGAR_SHIFT;
                    break;
                case BAM_CINS:
                    break;
                default:
                    j = p->b->core.n_cigar; // terminate for loop
                }
            }

            for (j = 0; j < 10 && ins_len[j]; j++)
                if (ins_len[j] == indel_len)
                    break;
            if (j < 10) {
                ins_len_freq[j]++;
                ins_len[j] = indel_len;
            }
            while (j > 0 && ins_len_freq[j] > ins_len_freq[j-1]) {
                int tmp = ins_len[j-1];
                ins_len[j-1] = ins_len[j];
                ins_len[j] = tmp;
                tmp = ins_len_freq[j-1];
                ins_len_freq[j-1] = ins_len_freq[j];
                ins_len_freq[j] = tmp;
            }
        }

        // Map ambiguity codes to one or more component bases.
        if (b < 16) {
            int Q = seqi2A[b] * (opts->use_qual ? q : 1);
            freq[1]  += Q?1:0;
            score[1] += Q?Q:0;
            Q = seqi2C[b] * (opts->use_qual ? q : 1);
            freq[2]  += Q?1:0;
            score[2] += Q?Q:0;
            Q = seqi2G[b] * (opts->use_qual ? q : 1);
            freq[4]  += Q?1:0;
            score[4] += Q?Q:0;
            Q = seqi2T[b] * (opts->use_qual ? q : 1);
            freq[8]  += Q?1:0;
            score[8] += Q?Q:0;
        } else { /* * */
            freq[16] ++;
            score[16]+=8 * (opts->use_qual ? q : 1);
        }
    }

    if (ins_len_freq[0] && nins >= nplp/2 && nins >= opts->min_depth) {
        if (insertion_consensus(plp, nplp, opts, nins, ins_len[0],
                                ins_ks, inq_ks) < 0)
            return -1;
    }

    // Total usable depth
    int tscore = 0;
    for (i = 0; i < 5; i++)
        tscore += score[1<<i];

    // Best and second best potential calls
    int call1  = 15, call2 = 15;
    int depth1 = 0,  depth2 = 0;
    int score1 = 0,  score2 = 0;
    for (i = 0; i < 5; i++) {
        int c = 1<<i; // A C G T *
        if (score1 < score[c]) {
            depth2 = depth1;
            score2 = score1;
            call2  = call1;
            depth1 = freq[c];
            score1 = score[c];
            call1  = c;
        } else if (score2 < score[c]) {
            depth2 = freq[c];
            score2 = score[c];
            call2  = c;
        }
    }

    // Work out which best and second best are usable as a call
    int used_score = score1;
    int used_depth = depth1;
    int used_base  = call1;
    if (score2 >= opts->het_fract * score1 && opts->ambig) {
        used_base  |= call2;
        used_score += score1;
        used_depth += depth2;
    }

    // N is too shallow, or insufficient proportion of total
    if (used_depth < opts->min_depth ||
        used_score < opts->call_fract * tscore) {
        used_depth = 0;
        // But note shallow gaps are still called gaps, not N, as
        // we're still more confident there is no base than it is
        // A, C, G or T.
        used_base = call1 == 16 /*&& depth1 >= call_fract * depth*/
            ? 16 : 0; // * or N
    }

    // Our final call.  "?" shouldn't be possible to generate
    const char *het =
        "NACMGRSVTWYHKDBN"
        "*ac?g???t???????";

    //printf("%c %d\n", het[used_base], used_depth);
    if (qual)
        *qual = 100.0 * used_score / tscore;
    return het[used_base];
}


/* --------------------------------------------------------------------------
 * Main processing logic
 */

// FIXME: move to header file if we intend to keep this interaction.
extern int pileup_seq(FILE *fp, const bam_pileup1_t *p, hts_pos_t pos,
                      hts_pos_t ref_len, const char *ref, kstring_t *ks,
                      int rev_del, int no_ins, int no_ins_mods,
                      int no_del, int no_ends);

int consensus_pileup(consensus_opts *opts, const bam_pileup1_t *p,
                     int np, int tid, int pos, int last_pos,
                     kstring_t *seq, kstring_t *qual) {
    kstring_t mod_ks = {0, 0}, ins_ks = {0, 0}, inq_ks = {0, 0};
    int cq, cb, cm_top = 0, cm_bot = 0, qm_top = 0, qm_bot = 0;
    char mod[10];

    int is_het = 0;
    if (opts->gap5) {
        consensus_t cons;
        if (opts->use_mqual)
            calculate_consensus_gap5(pos, CONS_MQUAL, p, np, opts,
                                     &cons, opts->default_qual,
                                     opts->mods ? &mod_ks : NULL,
                                     opts->mod_prob, opts->mod_yes,
                                     opts->mod_no, &ins_ks, &inq_ks);
        else
            calculate_consensus_gap5(pos, 0, p, np, opts,
                                     &cons, opts->default_qual,
                                     opts->mods ? &mod_ks : NULL,
                                     opts->mod_prob, opts->mod_yes,
                                     opts->mod_no, &ins_ks, &inq_ks);
        *mod = 0;
        if (cons.het_logodd > 0 && opts->ambig) {
            cb = "AMRWa" // 5x5 matrix with ACGT* per row / col
                 "MCSYc"
                 "RSGKg"
                 "WYKTt"
                 "acgt*"[cons.het_call];
            cq = cons.het_logodd;
            is_het = 1;
        } else{
            cb = "ACGT*"[cons.call];
            cq = cons.phred;

            // FIXME: no support for ChEBI here yet.
            if (cons.mod_qual[0] > opts->mod_cutoff) {
                cm_top = cons.mod_call[0];
                qm_top = cons.mod_qual[0];
                sprintf(mod,"+%c%d", cm_top, qm_top);
            }
            if (cons.mod_qual[1] > opts->mod_cutoff) {
                cm_bot = cons.mod_call[1];
                qm_bot = cons.mod_qual[1];
                sprintf(mod + strlen(mod),"-%c%d", cm_bot, qm_bot);
            }
        }
        if (cq < opts->cons_cutoff) {
            cb = 'N';
            cq = 0;
        }
    } else {
        cb = calculate_consensus_simple(p, np, opts, &cq, &ins_ks, &inq_ks);
        if (cb < 0)
            return -1;
    }

    if (seq) {
        if (pos > last_pos+1) {
            if (ks_expand(seq,  pos - last_pos) < 0 ||
                ks_expand(qual, pos - last_pos) < 0)
                return -1;
            memset(seq->s  + seq->l,  'N', pos - (last_pos+1));
            memset(qual->s + qual->l, '!', pos - (last_pos+1));
            seq->l  += pos - (last_pos+1);
            qual->l += pos - (last_pos+1);
        }

        // FIXME: no support for ChEBI here yet.
        // We need to permit a user-defined table of codes to chars.
        if (opts->fmt != PILEUP) {
            if (cm_top) {
                cb = cm_top;
                cq = qm_top;
            } else if (cm_bot) {
                switch (cm_bot) {
                case 'm': cb = '1'; break;
                case 'h': cb = '2'; break;
                case 'f': cb = '3'; break;
                case 'c': cb = '4'; break;
                // no complementary symbol avail.  Eg 'a'.  Replace with 'n'?
                default: cb = '?'; break;
                }
                cq = qm_bot;
            }
        }
        if (cb != '*' || opts->show_del) {
            kputc(cb, seq);

            kputc(MIN(cq, '~'-'!')+'!', qual);
        }
        if (ins_ks.l && opts->show_ins) {
            kputs(ins_ks.s, seq);
            kputs(inq_ks.s, qual);
        }
    }

    if (opts->fmt == PILEUP) {
        kstring_t ks = {0,0};
        int j;

        if (opts->het_only && (!is_het || islower(cb)))
            return 0;
        printf("%s\t%d\t%c%s\t%d\t",
               sam_hdr_tid2name(opts->h, tid), pos+1, cb, mod, cq);

        // FIXME: pileup_seq uses the mod iterator, but we've
        // already used it in consensus call.
        // Thus this doesn't report them as we're beyond that location.
        // FIXME: does this also break the next cons call?
        for (j = 0; j < np; j++)
            pileup_seq(stdout, p+j, pos, 0, NULL, &ks, 0, 2, 1, 2, 1);
            //pileup_seq(stdout, p+j, pos, 0, NULL, &ks, 0, 0, 0, 0, 0);
        putchar('\n');

        free(ks.s); // FIXME: reuse this
    }

    return 0;
}

static void dump_fastq(const char *name,
                       const char *seq, size_t seq_l,
                       const char *qual, size_t qual_l,
                       enum format fmt,
                       int line_len) {
    printf("%c%s\n", ">@"[fmt==FASTQ], name);
    size_t i;
    for (i = 0; i < seq_l; i += line_len)
        printf("%.*s\n", (int)MIN(line_len, seq_l - i), seq+i);

    if (fmt == FASTQ) {
        printf("+\n");
        for (i = 0; i < seq_l; i += line_len)
            printf("%.*s\n", (int)MIN(line_len, seq_l - i), qual+i);
    }
}

// Iterate over the pileup
static int consensus_loop(consensus_opts *opts) {
    bam_plp_t iter;
    const bam_pileup1_t *p;
    int tid, pos, n, last_tid = -99;
    kstring_t seq = {0}, qual = {0};
    int last_pos = -1;

    iter = bam_plp_init(readaln, (void *)opts);
    if (opts->mods) {
        bam_plp_constructor(iter, pileup_cd_create);
        bam_plp_destructor(iter, pileup_cd_destroy);
    }

    while ((p = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {
        if (opts->iter && (pos < opts->iter->beg || pos >= opts->iter->end))
            continue;
        if (tid != last_tid) {
            if (last_tid >= 0 && opts->fmt != PILEUP) {
                if (opts->all_bases) {
                    int i;
                    int N = opts->iter ? opts->iter->end : INT_MAX;
                    N = MIN(N, sam_hdr_tid2len(opts->h,last_tid)) - last_pos-1;
                    if (N > 0) {
                        ks_expand(&seq, N+1);
                        for (i = 0; i < N; i++)
                            seq.s[seq.l++] = 'N';
                        seq.s[seq.l] = 0;
                    }
                }
                dump_fastq(sam_hdr_tid2name(opts->h, last_tid),
                           seq.s, seq.l, qual.s, qual.l, opts->fmt,
                           opts->line_len);
            }
            seq.l = 0; qual.l = 0;
            last_tid = tid;
            if (opts->iter)
                last_pos = opts->iter->beg-1;
            else
                last_pos = opts->all_bases ? -1 : pos-1;
        }
        if (consensus_pileup(opts, p, n, tid, pos, last_pos, &seq, &qual) < 0)
            return -1;
        last_pos = pos;
    }
    bam_plp_destroy(iter);

    if (last_tid >= 0 && opts->fmt != PILEUP) {
        if (opts->all_bases) {
            int i;
            int N = opts->iter ? opts->iter->end : INT_MAX;
            N = MIN(N, sam_hdr_tid2len(opts->h, last_tid)) - last_pos - 1;
            if (N > 0) {
                ks_expand(&seq, N+1);
                for (i = 0; i < N; i++)
                    seq.s[seq.l++] = 'N';
                seq.s[seq.l] = 0;
            }
        }
        dump_fastq(sam_hdr_tid2name(opts->h, last_tid),
                   seq.s, seq.l, qual.s, qual.l, opts->fmt,
                   opts->line_len);
    }

    return 0;
}

static void usage_exit(FILE *fp, int exit_status) {
    fprintf(fp, "Usage: samtools consensus [options] <in.bam>\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "  -r, --region REG    Limit query to REG. Requires an index\n");
    fprintf(fp, "  -f, --format FMT    Output in format FASTA, FASTQ or PILEUP [FASTA]\n");
    fprintf(fp, "  -l, --line-len N    Wrap FASTA/Q at line length N [70]\n");
    fprintf(fp, "  -5                  Enable the bayesian 'gap5' consensus mode [off]\n");
    fprintf(fp, "  -a                  Output all bases (start/end of reference)\n");
    fprintf(fp, "  --show-del yes/no   Whether to show deletion as \"*\" [no]\n");
    fprintf(fp, "  --show-ins yes/no   Whether to show insertions [yes]\n");
    fprintf(fp, "For simple consensus mode:\n");
    fprintf(fp, "  -q, --use-qual      Use quality values in calculation\n");
    fprintf(fp, "  -d, --min-depth D   Minimum depth of D [1]\n");
    fprintf(fp, "  -c, --call-fract C  At least C portion of bases must agree [0.75]\n");
    fprintf(fp, "  -H, --het-fract C   Minimum fraction of 2nd-most to most common base [0.66]\n");
    fprintf(fp, "For gap5 consensus mode:\n");
    fprintf(fp, "  -C, --cutoff C      Consensus cutoff quality C [20]\n");
    fprintf(fp, "  -A, --ambig         Enable IUPAC ambiguity codes [off]\n");
    fprintf(fp, "  -M, --output-mods   Also report base modifications [off]\n");
    fprintf(fp, "  --mod-yes-cutoff x  P>=x to consider mod as true [0.8]\n");
    fprintf(fp, "  --mod-no-cutoff x   P<=x to consider mod as false [0.2]\n");
    fprintf(fp, "  --mod-fixed-prob p  Use probality P for all mods above cutoff [0.8]\n");
    fprintf(fp, "  --mod-cutoff Q      Emit mods if phred score >= Q [10]\n");

    // Plus gap5 vs simple
    // Plus -F format (fasta/fastq/pileup)
    fprintf(fp, "Global options:\n");
    sam_global_opt_help(fp, "-.---@-.");
    exit(exit_status);
}

int main_consensus(int argc, char **argv) {
    int c;

    consensus_opts opts = {
        .use_qual   = 0,
        .use_mqual  = 0,
        .min_depth  = 1,
        .call_fract = 0.75,
        .het_fract  = 0.66,
        .het_only   = 0,
        .fmt        = FASTA,
        .cons_cutoff= 20,
        .ambig      = 0,
        .line_len   = 70,
        .mods       = 0,
        .default_qual = 10,
        .mod_cutoff   = 10,
        .mod_prob     = 0.8,
        .mod_yes      = 0.8,
        .mod_no       = 0.2,
        .all_bases    = 0,
        .show_del     = 0,
        .show_ins     = 1,
    };

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', '-', '-', '@'),
        {"use-qual",           no_argument,       NULL, 'q'},
        {"use-mqual",          no_argument,       NULL, 'm'},
        {"min-depth",          required_argument, NULL, 'd'},
        {"call-fract",         required_argument, NULL, 'c'},
        {"het-fract",          required_argument, NULL, 'H'},
        {"region",             required_argument, NULL, 'r'},
        {"format",             required_argument, NULL, 'f'},
        {"cutoff",             required_argument, NULL, 'C'},
        {"ambig",              no_argument,       NULL, 'A'},
        {"line-len",           required_argument, NULL, 'l'},
        {"output-mods",        no_argument,       NULL, 'M'},
        {"default-qual",       required_argument, NULL, 1},
        {"mod-yes-cutoff",     required_argument, NULL, 2},
        {"mod-no-cutoff",      required_argument, NULL, 3},
        {"mod-fixed-prob",     required_argument, NULL, 4},
        {"mod-cutoff",         required_argument, NULL, 5},
        {"het-only",           no_argument,       NULL, 6},
        {"show-del",           required_argument, NULL, 7},
        {"show-ins",           required_argument, NULL, 8},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "@:qmd:c:H:r:5f:C:aAl:M",
                            lopts, NULL)) >= 0) {
        switch (c) {
        case 'a': opts.all_bases++; break;
        case 'q': opts.use_qual=1; break;
        case 'm': opts.use_mqual=1; break;
        case 'd': opts.min_depth = atoi(optarg); break;
        case 'c': opts.call_fract = atof(optarg); break;
        case 'H': opts.het_fract = atof(optarg); break;
        case 'r': opts.reg = optarg; break;
        case '5': opts.gap5 = 1; break;
        case 'C': opts.cons_cutoff = atoi(optarg); break;
        case 'A': opts.ambig = 1; break;
        case 'M': opts.mods = 1; break;
        case 1:   opts.default_qual = atoi(optarg); break;
        case 2:   opts.mod_yes = atof(optarg); break;
        case 3:   opts.mod_no  = atof(optarg); break;
        case 4:   opts.mod_prob = atof(optarg); break;
        case 5:   opts.mod_cutoff = atoi(optarg); break;
        case 6:   opts.het_only = 1; break;
        case 7:   opts.show_del = (*optarg == 'y' || *optarg == 'Y'); break;
        case 8:   opts.show_ins = (*optarg == 'y' || *optarg == 'Y'); break;
        case 'l':
            if ((opts.line_len = atoi(optarg)) <= 0)
                opts.line_len = INT_MAX;
            break;

        case 'f':
            if (strcasecmp(optarg, "fasta") == 0)
                opts.fmt = FASTA;
            else if (strcasecmp(optarg, "fastq") == 0)
                opts.fmt = FASTQ;
            else if (strcasecmp(optarg, "pileup") == 0)
                opts.fmt = PILEUP;
            else {
                fprintf(stderr, "Unknown format %s\n", optarg);
                return 1;
            }
            break;

        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    if (argc != optind+1) {
        if (argc == optind) usage_exit(stdout, EXIT_SUCCESS);
        else usage_exit(stderr, EXIT_FAILURE);
    }
    opts.fp = sam_open_format(argv[optind], "r", &ga.in);
    if (opts.fp == NULL) {
        print_error_errno("consensus", "Cannot open input file \"%s\"",
                          argv[optind]);
        return 1;
    }
    if (ga.nthreads > 0)
        hts_set_threads(opts.fp, ga.nthreads);

    if (hts_set_opt(opts.fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        return 1;
    }

    if (!(opts.h = sam_hdr_read(opts.fp))) {
        fprintf(stderr, "Failed to read header for \"%s\"\n", argv[optind]);
        return 1;
    }

    if (opts.reg) {
        opts.idx = sam_index_load(opts.fp, argv[optind]);
        if (!opts.idx) {
            print_error("consensus", "Cannot load index for input file \"%s\"",
                        argv[optind]);
            return 1;
        }
        opts.iter = sam_itr_querys(opts.idx, opts.h, opts.reg);
        if (!opts.iter) {
            print_error("consensus", "Failed to parse region \"%s\"",
                        opts.reg);
            return 1;
        }
    }

    if (consensus_loop(&opts) < 0) {
        print_error_errno("consensus", "Failed");
        return 1;
    }

    if (opts.iter)
        hts_itr_destroy(opts.iter);
    if (opts.idx)
        hts_idx_destroy(opts.idx);

    if (sam_close(opts.fp) < 0) {
        print_error_errno("consensus", "Closing input file \"%s\"",
                          argv[optind]);
        return 1;
    }

    sam_hdr_destroy(opts.h);
    sam_global_args_free(&ga);

    return 0;
}
