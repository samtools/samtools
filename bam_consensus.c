/*  bam_consensus.c -- consensus subcommand.

    Copyright (C) 1998-2001,2003 Medical Research Council (Gap4/5 source)
    Copyright (C) 2003-2005,2007-2023 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

The primary work here is GRL since 2021, under an MIT license.
Sections derived from Gap5, which include calculate_consensus_gap5()
associated functions, are mostly copyright Genome Research Limited from
2003 onwards.  These were originally under a BSD license, but as GRL is
copyright holder these portions can be considered to also be under the
same MIT license below:


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
 * The Gap5 consensus algorithm was in turn derived from the earlier Gap4
 * tool, developed by the Medical Research Council as part of the
 * Staden Package.  It is unsure how much of this source code is still
 * extant, without deep review, but the license used was a compatible
 * modified BSD license, included below.
 */

/*
Modified BSD license for any legacy components from the Staden Package:

Copyright (c) 2003 MEDICAL RESEARCH COUNCIL
All rights reserved

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   . Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

   . Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

   . Neither the name of the MEDICAL RESEARCH COUNCIL, THE LABORATORY OF
MOLECULAR BIOLOGY nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


// FIXME: also use strand to spot possible basecalling errors.
//        Specifically het calls where mods are predominantly on one
//        strand.  So maybe require + and - calls and check concordance
//        before calling a het as confident.  (Still call, but low qual?)

// TODO: call by kmers rather than individual bases?  Or use kmers to skew
// quality at least.  It can identify variants that are low quality due to
// neighbouring edits that aren't consistently correlated.

// TODO: pileup callback ought to know when it's the last in the region /
// chromosome.  This means the caller code doesn't have to handle the
// termination phase and deduplicates the code.  (Changing from
// one chr to the next is the same as ending the last.)
//
// TODO: track which reads contribute to multiple confirmed (HQ) differences
// vs which contribute to only one (LQ) difference.  Correlated changes
// are more likely to be real.  Ie consensus more of a path than solely
// isolated columns.
//
// Either that or a dummy "end of data" call is made to signify end to
// permit tidying up.  Maybe add a "start of data" call too?

// Eg 50T 20A seems T/A het,
// but 30T+ 20T- 18A+ 2A- seems like a consistent A miscall on one strand
// only, while T is spread evenly across both strands.

// TODO:  Phasing of long reads.
// Long reads offer very strong phasing opportunities for SNPs.
// From these, we get strong evidence for accuracy of indels.
// Specifically whether the distribution of poly-len within a phases
// is significantly different to the distribution of poly len between
// phases.

// TODO end STR trimming. Eg:
// REF AAGCTGAAAAGTTAATGTCTTATTTTTTTTTTTTTTTTGAGATGGAGTC
//     aagctgaaaagttaatgtctta****ttttttttttttgagatggagtc
//     aagctgaaaagttaatgtcttattttttttt
//     aagctgaaaagttaatgtctta****ttttttttttttgagatggagtc
// Middle seq doesn't validate those initial T alignments.
// Qual_train solves this by use of the STR trimmer.

// TODO add a weight for proximity to homopolymer.
// Maybe length/distance?  So 3 away from a 12-mer is similar to 1 away
// from a 4-mer?

// TODO: Count number of base types between this point and the nearest
// indel or end of read.  Eg GATCG<here>AGAGAG*TAGC => 2 (A and G).
// adj is nbase/4 * score, or (nbase+1)/5?
// Perhaps multiplied by length too, to get local complexity score?

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>

#include <htslib/sam.h>
#include <htslib/hfile.h>

#include "samtools.h"
#include "sam_opts.h"
#include "bam_plbuf.h"
#include "consensus_pileup.h"

#ifdef __SSE__
#   include <xmmintrin.h>
#else
#   define _mm_prefetch(a,b)
#endif

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef MAX
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

// Defines for experiment code which is currently disabled

// Hardy-Weinberg statistics to check heterozygous sites match allelic
// frequencies.
//#define DO_HDW

// Filter bayesian calls by min-depth and min-fract parameters
//#define DO_FRACT

// Checks uniqueness of surrounding bases to adjust scores
//#define K2 2

// Look for strand bias in distribution of homopolymer lengths
//#define DO_POLY_DIST

// Minimum cutoff for storing mod data; => at least 10% chance
#define MOD_CUTOFF 0.46

enum format {
    FASTQ,
    FASTA,
    PILEUP
};

typedef unsigned char uc;

// Simple recalibration table for substitutions, undercalls and overcalls.
// In future, we'll update this to be kmer based too.
typedef struct {
    int smap[101]; // substituion or SNP
    int umap[101]; // undercall   or DEL
    int omap[101]; // overcall    or INS
} qcal_t;

typedef struct {
    // User options
    char *reg;
    int use_qual;
    int min_qual;
    int adj_qual;
    int use_mqual;
    double scale_mqual;
    int nm_adjust;
    int nm_halo;
    int sc_cost;
    int low_mqual;
    int high_mqual;
    int min_depth;
    double call_fract;
    double het_fract;
    int mode;   // One of MODE_* macros below
    enum format fmt;
    int cons_cutoff;
    int ambig;
    int line_len;
    int default_qual;
    int het_only;
    int all_bases;
    int show_del;
    int show_ins;
    int mark_ins;
    int excl_flags;
    int incl_flags;
    int min_mqual;
    double P_het;
    double P_indel;
    double het_scale;
    double homopoly_fix;
    double homopoly_redux;
    qcal_t qcal;

    // Internal state
    samFile *fp;
    FILE *fp_out;
    sam_hdr_t *h;
    hts_idx_t *idx;
    hts_itr_t *iter;
    kstring_t ks_line;
    kstring_t ks_ins_seq;
    kstring_t ks_ins_qual;
    int last_tid;
    hts_pos_t last_pos;
} consensus_opts;

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
} consensus_t;

#define P_HET 1e-3
#define P_INDEL 2e-4
#define P_HOMOPOLY 0.5
#define P_HET_SCALE 1.0

#define LOG10            2.30258509299404568401
#define TENOVERLOG10     4.34294481903251827652
#define TENLOG2OVERLOG10 3.0103

#ifdef __GNUC__
#define ALIGNED(x) __attribute((aligned(x)))
#else
#define ALIGNED(x)
#endif

// Initialised once as a global array.  This won't work if threaded,
// but we'll rewrite if and when that gets added later.
static double e_tab_a[1002]  ALIGNED(16);
static double *e_tab = &e_tab_a[500];
static double e_tab2_a[1002] ALIGNED(16);
static double *e_tab2 = &e_tab2_a[500];
static double e_log[501]     ALIGNED(16);

/* Precomputed matrices for the consensus algorithm */
typedef struct {
    double prior[25]    ALIGNED(16);  /* Sum to 1.0 */
    double lprior15[15] ALIGNED(16);  /* 15 combinations of {ACGT*} */

    double pMM[101] ALIGNED(16);
    double p__[101] ALIGNED(16);
    double p_M[101] ALIGNED(16);
    double po_[101] ALIGNED(16);
    double poM[101] ALIGNED(16);
    double poo[101] ALIGNED(16);
    double puu[101] ALIGNED(16);
    double pum[101] ALIGNED(16);
    double pmm[101] ALIGNED(16);

    // Multiplier on homopolymer length before reducing phred qual
    double poly_mul;
} cons_probs;

// Two sets of params; recall oriented (gap5) and precision (stf).
// We use the former unless MODE_MIXED is set (which is the default
// for bayesian consensus mode if P_indel is significant).
static cons_probs cons_prob_recall, cons_prob_precise;

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
 *
 * We have a generic indel probability, but it's a catch all for overcall,
 * undercall, alignment artifacts, homopolymer issues, etc.  So we can set
 * it considerably higher and just let the QUAL skew do the filtering for
 * us, albeit no longer well calibrated.
 */

// NB: Should _M be MM?
// Ie sample really is A/C het, and we observe C.  That should be a match,
// not half a match.

#define MODE_SIMPLE    0 // freq counting

#define MODE_BAYES_116 1 // Samtools 1.16 (no indel param)
#define MODE_RECALL    2 // so called as it's the params from Gap5
#define MODE_PRECISE   3 // a more precise set; +FN, --FP
#define MODE_MIXED     4 // Combination of GAP5/BAYES

#define QCAL_FLAT           0
#define QCAL_HIFI           1
#define QCAL_HISEQ          2
#define QCAL_ONT_R10_4_SUP  3
#define QCAL_ONT_R10_4_DUP  4
#define QCAL_ULTIMA         5

// Calibration tables here don't necessarily reflect the true accuracy.
// They have been manually tuned to work in conjunction with other command
// line parameters used in the machine profiles.  For example reducing one
// qual here and increasing sensitivity elsewhere via another parameter.
static qcal_t static_qcal[6] = {
    { // FLAT
        {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
         10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
         20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
         30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
         40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
         50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
         60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
         70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
         80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
         90, 91, 92, 93, 94, 95, 96, 97, 98, 99},
        {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
         10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
         20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
         30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
         40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
         50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
         60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
         70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
         80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
         90, 91, 92, 93, 94, 95, 96, 97, 98, 99},
        {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
         10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
         20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
         30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
         40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
         50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
         60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
         70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
         80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
         90, 91, 92, 93, 94, 95, 96, 97, 98, 99}
    },

    { // HiFi
        {10, 11, 11, 12, 13, 14, 15, 16, 18, 19,
         20, 21, 22, 23, 24, 25, 27, 28, 29, 30,
         31, 32, 33, 33, 34, 35, 36, 36, 37, 38,
         38, 39, 39, 40, 40, 41, 41, 41, 41, 42,
         42, 42, 42, 43, 43, 43, 43, 43, 43, 43,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         },
        {  4,  4,  4,  4,  5,  6,  6,  7,  8,  9,
          10, 11, 11, 12, 13, 14, 15, 15, 16, 17,
          18, 19, 19, 20, 20, 21, 22, 23, 23, 24,
          25, 25, 25, 26, 26, 26, 27, 27, 28, 28,
          28, 28, 27, 27, 27, 28, 28, 28, 28, 27,
          27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
          27, 27, 26, 26, 25, 26, 26, 27, 27, 27,
          26, 26, 26, 26, 26, 26, 26, 26, 27, 27,
          28, 29, 28, 28, 28, 27, 27, 27, 27, 27,
          27, 28, 28, 30, 30, 30, 30, 30, 30, 30,
          },
        {  8,  8,  8,  8,  9, 10, 11, 12, 13, 14,
          15, 15, 16, 17, 18, 19, 19, 20, 20, 21,
          21, 22, 22, 23, 23, 23, 24, 24, 24, 25,
          25, 25, 25, 25, 25, 26, 26, 26, 26, 27,
          27, 27, 27, 27, 27, 28, 28, 28, 28, 28,
          29, 29, 29, 29, 29, 29, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
        }
    },

    { // HiSeq
        { 2,  2,  2,  3,  3,  4,  5,  5,  6,  7,
          8,  9, 10, 11, 11, 12, 13, 14, 15, 16,
          17, 17, 18, 19, 20, 21, 22, 22, 23, 24,
          25, 26, 27, 28, 28, 29, 30, 31, 32, 33,
          34, 34, 35, 36, 37, 38, 39, 39, 40, 41,
          42, 43, 44, 45, 45, 46, 47, 48, 49, 50,
          51, 51, 52, 53, 54, 55, 56, 56, 57, 58,
          59, 60, 61, 62, 62, 63, 64, 65, 66, 67,
          68, 68, 69, 70, 71, 72, 73, 73, 74, 75,
          76, 77, 78, 79, 79, 80, 81, 82, 83, 84,
        },
        { 1,  2,  3,  4,  5,  7,  8,  9, 10, 11,
          13, 14, 15, 16, 17, 19, 20, 21, 22, 23,
          25, 26, 27, 28, 29, 31, 32, 33, 34, 35,
          37, 38, 39, 40, 41, 43, 44, 45, 46, 47,
          49, 50, 51, 52, 53, 55, 56, 57, 58, 59,
          61, 62, 63, 64, 65, 67, 68, 69, 70, 71,
          73, 74, 75, 76, 77, 79, 80, 81, 82, 83,
          85, 86, 87, 88, 89, 91, 92, 93, 94, 95,
          97, 98, 99, 100, 101, 103, 104, 105, 106, 107,
          109, 110, 111, 112, 113, 115, 116, 117, 118, 119,
        },
        { 1,  2,  3,  4,  5,  7,  8,  9, 10, 11,
          13, 14, 15, 16, 17, 19, 20, 21, 22, 23,
          25, 26, 27, 28, 29, 31, 32, 33, 34, 35,
          37, 38, 39, 40, 41, 43, 44, 45, 46, 47,
          49, 50, 51, 52, 53, 55, 56, 57, 58, 59,
          61, 62, 63, 64, 65, 67, 68, 69, 70, 71,
          73, 74, 75, 76, 77, 79, 80, 81, 82, 83,
          85, 86, 87, 88, 89, 91, 92, 93, 94, 95,
          97, 98, 99, 100, 101, 103, 104, 105, 106, 107,
          109, 110, 111, 112, 113, 115, 116, 117, 118, 119,
        }
    },
    { // ONT R10.4 super
        {  0,  2,  2,  2,  3,  4,  4,  5,  6,  7,
           7,  8,  9, 12, 13, 14, 15, 15, 16, 17,
          18, 19, 20, 22, 24, 25, 26, 27, 28, 29,
          30, 31, 33, 34, 36, 37, 38, 38, 39, 39,
          40, 40, 40, 40, 40, 40, 40, 41, 40, 40,
          41, 41, 40, 40, 40, 40, 41, 40, 40, 40,
          40, 41, 41, 40, 40, 41, 40, 40, 39, 41,
          40, 41, 40, 40, 41, 41, 41, 40, 40, 40,
          40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
          40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        },
        {  0,  2,  2,  2,  3,  4,  5,  6,  7,  8,
           8,  9,  9, 10, 10, 10, 11, 12, 12, 13,
          13, 13, 14, 14, 15, 16, 16, 17, 18, 18,
          19, 19, 20, 21, 22, 23, 24, 25, 25, 25,
          25, 25, 25, 25, 25, 25, 26, 26, 26, 26,
          26, 26, 26, 26, 27, 27, 27, 27, 27, 27,
          27, 27, 27, 27, 27, 27, 27, 28, 28, 28,
          28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
          28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
          28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
        },
        {  0,  4,  6,  6,  6,  7,  7,  8,  9,  9,
           9, 10, 10, 11, 11, 12, 12, 13, 13, 14,
          15, 15, 15, 16, 16, 17, 17, 18, 18, 19,
          19, 20, 20, 21, 22, 22, 23, 23, 24, 24,
          24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
          24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
          24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
          24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
          24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
          24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
        }
    },
    { // ONT R10.4 duplex; just a copy of hifi for now
        {10, 11, 11, 12, 13, 14, 15, 16, 18, 19,
         20, 21, 22, 23, 24, 25, 27, 28, 29, 30,
         31, 32, 33, 33, 34, 35, 36, 36, 37, 38,
         38, 39, 39, 40, 40, 41, 41, 41, 41, 42,
         42, 42, 42, 43, 43, 43, 43, 43, 43, 43,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
         },
        {  4,  4,  4,  4,  5,  6,  6,  7,  8,  9,
          10, 11, 11, 12, 13, 14, 15, 15, 16, 17,
          18, 19, 19, 20, 20, 21, 22, 23, 23, 24,
          25, 25, 25, 26, 26, 26, 27, 27, 28, 28,
          28, 28, 27, 27, 27, 28, 28, 28, 28, 27,
          27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
          27, 27, 26, 26, 25, 26, 26, 27, 27, 27,
          26, 26, 26, 26, 26, 26, 26, 26, 27, 27,
          28, 29, 28, 28, 28, 27, 27, 27, 27, 27,
          27, 28, 28, 30, 30, 30, 30, 30, 30, 30,
          },
        {  8,  8,  8,  8,  9, 10, 11, 12, 13, 14,
          15, 15, 16, 17, 18, 19, 19, 20, 20, 21,
          21, 22, 22, 23, 23, 23, 24, 24, 24, 25,
          25, 25, 25, 25, 25, 26, 26, 26, 26, 27,
          27, 27, 27, 27, 27, 28, 28, 28, 28, 28,
          29, 29, 29, 29, 29, 29, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
          30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
        }
    },
    { // Ultima Genomics
        { 2,  2,  3,  4,  5,  6,  6,  7,  8,  9,
          10, 10, 11, 12, 13, 14, 14, 15, 16, 17,
          18, 18, 19, 21, 22, 23, 23, 24, 25, 26,
          27, 27, 28, 29, 30, 31, 31, 32, 33, 34,
          35, 35, 36, 37, 38, 39, 39, 40, 42, 43,
          44, 44, 45, 46, 47, 48, 48, 49, 50, 51,
          52, 52, 53, 54, 55, 56, 56, 57, 58, 59,
          60, 60, 61, 63, 64, 65, 65, 66, 67, 68,
          69, 69, 70, 71, 72, 73, 73, 74, 75, 76,
          77, 77, 78, 79, 80, 81, 81, 82, 84, 85,
        },
        { 1,  1,  2,  2,  3,  3,  4,  4,  4,  4,
          5,  5,  6,  6,  7,  7,  8,  8,  9, 10,
          10, 10, 11, 12, 13, 13, 13, 14, 15, 16,
          16, 16, 17, 18, 18, 19, 19, 20, 20, 21,
          21, 22, 22, 22, 22, 23, 23, 24, 24, 25,
          25, 25, 25, 25, 25, 25, 26, 26, 26, 26,
          26, 26, 27, 27, 27, 27, 27, 27, 27, 27,
          27, 28, 28, 28, 28, 28, 28, 28, 28, 28,
          28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
          28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
        },
        { 1,  1,  2,  2,  3,  3,  4,  4,  4,  4,
          5,  5,  6,  6,  7,  7,  8,  8,  9, 10,
          10, 10, 11, 12, 13, 13, 13, 14, 15, 16,
          16, 16, 17, 18, 18, 19, 19, 20, 20, 21,
          21, 22, 22, 22, 22, 23, 23, 24, 24, 25,
          25, 25, 25, 25, 25, 25, 26, 26, 26, 26,
          26, 26, 27, 27, 27, 27, 27, 27, 27, 27,
          27, 28, 28, 28, 28, 28, 28, 28, 28, 28,
          28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
          28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
        }
    }
};

int set_qcal(qcal_t *q, int id) {
    if (id < 0 || id >= sizeof(static_qcal)/sizeof(*static_qcal))
        return -1;

    memcpy(q, &static_qcal[id], sizeof(*q));
    return 0;
}

int load_qcal(qcal_t *q, const char *fn) {
    int i;

    if (strcmp(fn, ":hifi") == 0)
        return set_qcal(q, QCAL_HIFI);
    if (strcmp(fn, ":hiseq") == 0)
        return set_qcal(q, QCAL_HISEQ);
    if (strcmp(fn, ":r10.4_sup") == 0)
        return set_qcal(q, QCAL_ONT_R10_4_SUP);
    if (strcmp(fn, ":r10.4_dup") == 0)
        return set_qcal(q, QCAL_ONT_R10_4_DUP);
    if (strcmp(fn, ":ultima") == 0)
        return set_qcal(q, QCAL_ULTIMA);

    // default
    for (i = 0; i < 101; i++)
        q->smap[i] = q->umap[i] = q->omap[i] = i;

    if (strcmp(fn, ":flat") == 0)
        return 0;

    hFILE *fp = hopen(fn, "r");
    if (!fp)
        return -1;

    kstring_t line = KS_INITIALIZE;
    int max = 0;
    int last_qual = 0;
    while (line.l = 0, kgetline(&line, (kgets_func *)hgets, fp) >= 0) {
        int v, s, u, o;
        if (*line.s == '#')
            continue;
        if (sscanf(line.s, "QUAL %d %d %d %d", &v, &s, &u, &o) != 4)
            goto err;
        while (v > last_qual) {
            q->smap[last_qual+1] = q->smap[last_qual];
            q->umap[last_qual+1] = q->umap[last_qual];
            q->omap[last_qual+1] = q->omap[last_qual];
            last_qual++;
        }
        if (v >= 0 && v < 100) {
            q->smap[v] = s;
            q->umap[v] = u;
            q->omap[v] = o;
        }
        if (v < max) {
            fprintf(stderr, "Qual calibration file is not in ascending order\n");
            return hclose(fp) ? -2 : -1;
        }
        max = v;
    }

    for (i = max+1; i < 101; i++) {
        q->smap[i] = q->smap[max];
        q->umap[i] = q->umap[max];
        q->omap[i] = q->omap[max];
    }

    ks_free(&line);
    return hclose(fp) < 0 ? -2 : 0;

 err:
    ks_free(&line);
    return hclose(fp) < 0 ? -2 : -1;
}

static void consensus_init(double p_het, double p_indel, double het_scale,
                           double poly_mul,
                           qcal_t *qcal, int mode, cons_probs *cp) {
    int i;

    // NB: only need to initialise once, but we do here for now
    for (i = -500; i <= 500; i++)
        e_tab[i] = exp(i);
    for (i = -500; i <= 500; i++)
        e_tab2[i] = exp(i/10.);
    for (i = 0; i <= 500; i++)
        e_log[i] = log(i);

    // EXPERIMENTAL
    cp->poly_mul = poly_mul;

    // The priors make very little difference, unless shallow data.
    // ACGT* by ACGT*
    // So AA=0, CC=6, GG=12, TT=18, **=24
    for (i = 0; i < 25; i++)
        cp->prior[i] = p_het / 6; // AC AG AT CG CT GT

    // Flat assumption that it is what we observe, and measure everything else
    // as relative to this.
    cp->prior[0]=cp->prior[6]=cp->prior[12]=cp->prior[18]=cp->prior[24] = 1;

    // heterozygous deletion
    for (i = 4; i < 24; i+=5)
        cp->prior[i] = p_indel / 6; // /6 to be scaled vs p_het equivalently

    // heterozygous insertion
    for (i = 20; i < 24; i++)
        cp->prior[i] = p_indel / 6;

    cp->lprior15[0]  = log(cp->prior[0]);
    cp->lprior15[1]  = log(cp->prior[1]);
    cp->lprior15[2]  = log(cp->prior[2]);
    cp->lprior15[3]  = log(cp->prior[3]);
    cp->lprior15[4]  = log(cp->prior[4]);
    cp->lprior15[5]  = log(cp->prior[6]);
    cp->lprior15[6]  = log(cp->prior[7]);
    cp->lprior15[7]  = log(cp->prior[8]);
    cp->lprior15[8]  = log(cp->prior[9]);
    cp->lprior15[9]  = log(cp->prior[12]);
    cp->lprior15[10] = log(cp->prior[13]);
    cp->lprior15[11] = log(cp->prior[14]);
    cp->lprior15[12] = log(cp->prior[18]);
    cp->lprior15[13] = log(cp->prior[19]);
    cp->lprior15[14] = log(cp->prior[24]);

    for (i = 1; i < 101; i++) {
        double prob = 1 - pow(10, -qcal->smap[i] / 10.0);

        // Or is it that prob is 1-p(subst)-p(overcall)?
        cp->pMM[i] = log(prob);

        //cp->p__[i] = log(1-prob); // Big help to PB-CCS SNPs; unless fudged
        cp->p__[i] = log((1-prob)/3); // correct? poor on PB-CCS w/o fudge

        // Mixed alleles; just average two likelihoods
        cp->p_M[i] = log((exp(cp->pMM[i]) + exp(cp->p__[i]))/2);

        // What does this really mean?  Can we simulate this by priors?
        // It reduces the likelihood of calling het sites, which is
        // maybe compensation for alignment artifacts?  I'm unsure,
        // but it works (to differing degrees) on both PacBio HiFi and
        // Illumina HiSeq.  It (obviously) loses true hets, but
        // potentially this can be compensated for by tweaking P-het
        // (which is entirely in the priors).
        //
        // Low het_scale reduces false positives by making hets less
        // likely to be called.  In  high depth data we normally  have
        // enough evidence to call correctly even with low het_scale,
        // so it's a good +FN vs --FP tradeoff.  However on low depth
        // data, het_scale can filter out too many true variants.
        //
        // TODO: So consider adjusting at the end maybe?
        // Also consider never changing calls, but changing their
        // confidence, so the data is what produces the call with the
        // parameters skewing the quality score distribution.
        cp->p_M[i] += log(het_scale);

        if (mode == MODE_BAYES_116) {
            // Compatibility with samtools 1.16

            // This had no differention for indel vs substitution error rates,
            // so o(vercall) and u(undercall) are subst(_).
            cp->pmm[i] = cp->pMM[i];
            cp->poM[i] = cp->p_M[i];
            cp->pum[i] = cp->p_M[i];
            cp->po_[i] = cp->p__[i];
            cp->poo[i] = cp->p__[i];
            cp->puu[i] = cp->p__[i];

        } else {
            // When observing A C G T; leads to insertion calls
            prob = 1 - pow(10, -qcal->omap[i] / 10.0);
            // /3 for consistency with ACGT rem as relative likelihoods.
            // Otherwise with flat priors we end up calling all shallow data
            // as "*", which is illogical.
            cp->poo[i] = log((1-prob)/3);

            // Ensure pMM is always more likely. (NB: This shouldn't happen
            // now with the addition of the /3 step above.)
            if (cp->poo[i] > cp->pMM[i]-.5)
                cp->poo[i] = cp->pMM[i]-.5;

            cp->po_[i] = log((exp(cp->poo[i]) + exp(cp->p__[i]))/2);
            cp->poM[i] = log((exp(cp->poo[i]) + exp(cp->pMM[i]))/2);

            // Overcalls should never be twice as likely than mismatches.
            // Het bases are mix of _M (other) and MM ops (this).
            // It's fine for _M to be less likely than oM (more likely
            // to be overcalled than miscalled),  but it should never
            // be stronger when combined with other mixed data.
            if (cp->poM[i] > cp->p_M[i]+.5)
                cp->poM[i] = cp->p_M[i]+.5;

            // Note --low-MQ and --scale-MQ have a big impact on
            // undercall errs.  May need to separate these options per
            // type, but how?
            // Multiple-calls, as with mixed mode?  This feels like a cheat

            prob = 1 - pow(10, -qcal->umap[i] / 10.0);
            cp->pmm[i] = log(prob);
            cp->puu[i] = log((1-prob)/3);
            if (cp->puu[i] > cp->pMM[i]-.5) // MM is -ve
                cp->puu[i] = cp->pMM[i]-.5;

            cp->pum[i] = log((exp(cp->puu[i]) + exp(cp->pmm[i]))/2);
        }
    }

    cp->pMM[0] = cp->pMM[1];
    cp->p__[0] = cp->p__[1];
    cp->p_M[0] = cp->p_M[1];

    cp->pmm[0] = cp->pmm[1];
    cp->poo[0] = cp->poo[1];
    cp->po_[0] = cp->po_[1];
    cp->poM[0] = cp->poM[1];
    cp->puu[0] = cp->puu[1];
    cp->pum[0] = cp->pum[1];
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

/* Taylor (deg 3) implementation of the log */
static inline double fast_log2(double val)
{
    // FP representation is exponent & mantissa, where
    // value = 2^E * M.
    // Hence log2(value) = log2(2^E * M)
    //                   = log2(2^E)+ log2(M)
    //                   =        E + log2(M)
    union { double d; uint64_t x; } u = {val};
    const int E = ((u.x >> 52) & 2047) - 1024; // exponent E
    // Initial log2(M) based on mantissa
    u.x &= ~(2047LL << 52);
    u.x +=   1023LL << 52;

    val = ((-1/3.) * u.d + 2) * u.d - 2/3.;

    return E + val;
}

#define ph_log(x) (-TENLOG2OVERLOG10*fast_log2((x)))


int nins(const bam1_t *b){
    int i, indel = 0;
    uint32_t *cig = bam_get_cigar(b);
    for (i = 0; i < b->core.n_cigar; i++) {
        int op = bam_cigar_op(cig[i]);
        if (op == BAM_CINS || op == BAM_CDEL)
            indel += bam_cigar_oplen(cig[i]);
    }
    return indel;
}

/*
 * Some machines, including 454 and PacBio, store the quality values in
 * homopolymers with the first or last base always being the low quality
 * state.  This can cause problems when reverse-complementing and aligning,
 * especially when we left-justify indels.
 *
 * Other platforms take the approach of having the middle bases high and
 * the low confidence spread evenly to both start and end.  This means
 * reverse-complementing doesn't introduce any strand bias.
 *
 * We redistribute qualities within homopolymers in this style to fix
 * naive consensus or variant calling algorithms.
 */
void homopoly_qual_fix(bam1_t *b) {
    static double ph2err[256] = {0};
    int i;
    if (!ph2err[0]) {
        for (i = 0; i < 256; i++)
            ph2err[i] = pow(10, i/-10.0);
    }
    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);
    for (i = 0; i < b->core.l_qseq; i++) {
        int s = i; // start of homopoly
        int base = bam_seqi(seq, i);
        while (i+1 < b->core.l_qseq && bam_seqi(seq, i+1) == base)
            i++;
        // s..i inclusive is now homopolymer

        if (s == i)
            continue;

        // Simplest:  reverse if end_qual < start_qual
        // Next:      average outer-most two, then next two, etc
        // Best:      fully redistribute so start/end lower qual than centre

        // Middle route of averaging outer pairs is sufficient?
        int j, k;
        for (j = s, k = i; j < k; j++,k--) {
            double e = ph2err[qual[j]] + ph2err[qual[k]];
            qual[j] = qual[k] = -fast_log2(e/2)*3.0104+.49;
        }
    }
}

// Return the local NM figure within halo (+/- HALO) of pos.
// This local NM is used as a way to modify MAPQ to get a localised MAPQ
// score via an adhoc fashion.
double nm_local(const pileup_t *p, const bam1_t *b, hts_pos_t pos) {
    int *nm = (int *)p->cd;
    if (!nm)
        return 0;
    pos -= b->core.pos;
    if (pos < 0)
        return nm[0] & ((1<<24)-1);
    if (pos >= b->core.l_qseq)
        return nm[b->core.l_qseq-1] & ((1<<24)-1);

    return (nm[pos] & ((1<<24)-1)) / 10.0;
}

int poly_len(const pileup_t *p, const bam1_t *b, hts_pos_t pos) {
    int *nm = (int *)p->cd;
    if (!nm)
        return 0;
    pos -= b->core.pos;
    if (pos >= 0 && pos < b->core.l_qseq)
        return nm[pos] >> 24;
    else
        return 0;
}

/*
 * Initialise a new sequence appearing in the pileup.  We use this to
 * precompute some metrics that we'll repeatedly use in the consensus
 * caller; the localised NM score.
 *
 * We also directly amend the BAM record (which will be discarded later
 * anyway) to modify qualities to account for local quality minima.
 *
 * Returns 0 (discard) or 1 (keep) on success, -1 on failure.
 */
int nm_init(void *client_data, samFile *fp, sam_hdr_t *h, pileup_t *p) {
    consensus_opts *opts = (consensus_opts *)client_data;
    if (!opts->use_mqual)
        return 1;

    const bam1_t *b = &p->b;
    int qlen = b->core.l_qseq, i;
    if (qlen <= 0)
        return 0;
    int *local_nm = calloc(qlen, sizeof(*local_nm));
    if (!local_nm)
        return -1;
    p->cd = local_nm;

    double poly_adj = opts->homopoly_fix ? opts->homopoly_fix : 1;

    if (opts->adj_qual) {
        // Set local_nm based on a function of current qual and the local
        // minimum qual within the surrounding window.
        //
        // Basically if we're in a region of low confidence then we downgrade
        // higher qual outliers as they may not be as trustworthy as they
        // claim.  This may be because the qualities have been assigned to
        // the wrong or arbitrary base (very common in homopolymers), or the
        // surrounding quality (hence also error likelihood) have lead to
        // misalignments and the base may be contributing to the wrong
        // pileup column.
        //
        // The nm_local() function returns these scores and uses it to bias
        // the mapping quality, which in turn adjusts base quality.
        uint8_t *qual = bam_get_qual(b);
        uint8_t *seq = bam_get_seq(b);
        const int qhalo = 8; // window size for base qual
        int qmin = qual[0];  // min qual within qhalo
        const int qhalop = 2;// window size for homopolymer qual
        int qminp = qual[0]; // min qual within homopolymer halo
        int base = bam_seqi(seq, 0), polyl = 0, polyr = 0; // pos, not len

        // Minimum quality of the initial homopolymer
        for (i = 1; i < qlen; i++) {
            if (bam_seqi(seq, i) != base)
                break;
            if (i < qhalop && qminp > qual[i])
                qminp = qual[i];
        }

        // Minimum quality for general bases
        for (i = 0; i < qlen && i < qhalo; i++) {
            if (qmin > qual[i])
                qmin = qual[i];
        }

        for (;i < qlen-qhalo; i++) {
            if (opts->homopoly_fix && bam_seqi(seq, i) != base) {
                polyl = i;
                base = bam_seqi(seq, i);
                qminp = qual[i];
                int j;
                for (j = i+1; j < qlen; j++) {
                    if (bam_seqi(seq, j) != base)
                        break;
                    if (i < qhalop && qminp > qual[j])
                        qminp = qual[j];
                }
                polyr = j-1;
            } else {
                // CHECK: do we want to have opts->homopoly_fix above,
                // so when not in use we don't define pl to non-zero?
                // Test on SynDip
                polyr = polyl;
            }
            int pl = polyr-polyl;

            // Useful for SNPS, as we're judging the variation in
            // length as an indicator for base-misalignment.
            // Not so useful for indel calling where the longer the indel
            // the less confident we are on the len variation being real.
            int t = (opts->mode == MODE_BAYES_116)
                ? (qual[i]   + 5*qmin)/4
                : qual[i]/3 + (qminp-pl*2)*poly_adj;


            local_nm[i] += t < qual[i] ? qual[i]-t : 0;

            // Brute force qminp in polyl to polyr range.
            // TODO: optimise this with sliding window
            qminp = qual[i];
            int k;
            for (k = MAX(polyl,i-qhalop); k <= MIN(polyr,i+qhalop); k++)
                if (qminp > qual[k])
                    qminp = qual[k];

            if (qmin > qual[i+qhalo])
                qmin = qual[i+qhalo];
            else if (qmin <= qual[i-qhalo]) {
                int j;
                qmin = 99;
                for (j = i-qhalo+1; j <= i+qhalo; j++)
                    if (qmin > qual[j])
                        qmin = qual[j];
            }
        }
        for (; i < qlen; i++) {
            int t = (opts->mode == MODE_BAYES_116)
                ? (qual[i]   + 5*qmin)/4
                : qual[i]/3 + qminp*poly_adj;
            local_nm[i] += t < qual[i] ? qual[i]-t : 0;
        }
    }

    // Fix e.g. PacBio homopolymer qualities
    if (opts->homopoly_fix)
        homopoly_qual_fix((bam1_t *)b);

    // local_nm[i] & ((1<<24)-1) is for SNP score adjustment.
    // We also put some more basic poly-X len in local_nm[i] >> 24.
    uint8_t *seq = bam_get_seq(b);
    for (i = 0; i < qlen; i++) {
        int base = bam_seqi(seq, i);
        int poly = 0, j, k;
        for (j = i+1; j < qlen; j++)
            if (bam_seqi(seq, j) != base)
                break;
        //printf("%d x %d\n", base, j-i);

        poly = j-i-1; if (poly > 100) poly = 100;
        const int HALO=0;
        for (k = i-HALO; k < j+HALO; k++)
            if (k >= 0 && k < qlen)
                local_nm[k] = ((MAX(poly, local_nm[k]>>24))<<24)
                            | (local_nm[k] & ((1<<24)-1));

        i = j-1;
    }

    // Adjust local_nm array by the number of edits within
    // a defined region (pos +/- halo).
    const int halo = opts->nm_halo;
    const uint8_t *md = bam_aux_get(b, "MD");
    if (!md)
        return 1;
    md = (const uint8_t *)bam_aux2Z(md);

    // Handle cost of being near a soft-clip
    uint32_t *cig = bam_get_cigar(b);
    int ncig = b->core.n_cigar;

    if ( (cig[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP ||
        ((cig[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP && ncig > 1 &&
         (cig[1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)) {
        for (i = 0; i < halo && i < qlen; i++)
            local_nm[i]+=opts->sc_cost;
        for (; i < halo*2 && i < qlen; i++)
            local_nm[i]+=opts->sc_cost>>1;
    }
    if ( (cig[ncig-1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP ||
        ((cig[ncig-1] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP && ncig > 1 &&
         (cig[ncig-2] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)) {
        for (i = qlen-1; i >= qlen-halo && i >= 0; i--)
            local_nm[i]+=opts->sc_cost;
        for (; i >= qlen-halo*2 && i >= 0; i--)
            local_nm[i]+=opts->sc_cost>>1;
    }

    // Now iterate over MD tag
    int pos = 0;
    while (*md) {
        if (isdigit(*md)) {
            uint8_t *endptr;
            long i = strtol((char *)md, (char **)&endptr, 10);
            md = endptr;
            pos += i;
            continue;
        }

        // deletion.
        // Should we bump local_nm here too?  Maybe
        if (*md == '^') {
            while (*++md && !isdigit(*md))
                continue;
            continue;
        }

        // substitution
        for (i = pos-halo*2 >= 0 ?pos-halo*2 :0; i < pos-halo && i < qlen; i++)
            local_nm[i]+=5;
        for (; i < pos+halo && i < qlen; i++)
            local_nm[i]+=10;
        for (; i < pos+halo*2 && i < qlen; i++)
            local_nm[i]+=5;
        md++;
    }

    return 1;
}

void nm_free(void *client_data, samFile *fp, sam_hdr_t *h, pileup_t *p) {
    free(p->cd);
    p->cd = NULL;
}

#ifdef DO_HDW
/*
 * Stirling's formula with a 1/12n correction applied to improve accuracy.
 * This seems to hold remarkably true for both low and high numbers too.
 */
double lnfact(double n) {
    /* Or Gosper's formula...
     * return (n*ln(n) - n + ln(2*M_PI*n + M_PI/3) / 2);
     */
    return ((n+0.5)*log(n) - n + log(2*M_PI)/2) + log(1 + 1/(12.0*n));
        /* + log(1 + 1/(288.0*n*n)); */
}

/*
 * The binomical coefficient (n,k) for n trials with k successes where
 * prob(success) = p.
 *                               k      n-k
 * P (k|n) = n! / (k! (n-k)!)   p  (1-p)
 *  p
 *
 * The coefficient we are returning here is the n! / (k! (n-k)!) bit.
 * We compute it using ln(n!) and then exp() the result back to avoid
 * excessively large numbers.
 */
double bincoef(int n, double k) {
    return exp(lnfact(n) - lnfact(k) - lnfact(n-k));
}

/*
 * Given p == 0.5 the binomial expansion simplifies a bit, so we have
 * a dedicated function for this.
 */
double binprobhalf(int n, double k) {
    return bincoef(n, k) * pow(0.5, n);
}

double lnbinprobhalf(int n, double k) {
    // ln(binprobhalf) expanded up and simplified
    return lnfact(n) - lnfact(k) - lnfact(n-k) - 0.69315*n;
}
#endif

static
int calculate_consensus_gap5(hts_pos_t pos, int flags, int depth,
                             pileup_t *plp, consensus_opts *opts,
                             consensus_t *cons, int default_qual,
                             cons_probs *cp) {
    int i, j;
    static int init_done =0;
    static double q2p[101], mqual_pow[256];
    double min_e_exp = DBL_MIN_EXP * log(2) + 1;

    double S[15] ALIGNED(16) = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double sumsC[6] = {0,0,0,0,0,0}; // A C G T * N

    // Small hash on seq to check for uniqueness of surrounding bases.
    // If it's frequent, then it's more likely to be correctly called than
    // if it's rare.
    // Helps a bit on deep data, especially with K2=3, but detrimental on
    // shallow and (currently) quite a slow down.
#ifdef K2
    int hashN[1<<(K2*4+2)] = {0};
    int hash1[1<<2] = {0};
#endif

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
#ifdef DO_FRACT
    int counts2[2][6] = {{0}};
#endif

    /* Accumulate */

#ifdef K2
    const pileup_t *ptmp = plp;
    for (; ptmp; ptmp = ptmp->next) {
        const pileup_t *p = ptmp;
        if (p->qual < opts->min_qual)
            continue;

        int hb = 0;
#define _ 0
        static int X[16] = {_,0,1,_,2,_,_,_,3,_,_,_,_,_,_,_};
#undef _
        uint8_t *seq = bam_get_seq(&p->b);
        int i, base1 = X[p->base4];
        hash1[base1]++;
        for (i = p->seq_offset-K2; i <= p->seq_offset+K2; i++) {
            int base = i >= 0 && i < p->b.core.l_qseq ? X[bam_seqi(seq,i)] : _;
            hb = (hb<<2)|base;
        }
        hashN[hb]++;
    }
#endif

    int td = depth; // original depth
    depth = 0;
#ifdef DO_POLY_DIST
    int poly_dist[2][100] = {0};
#endif
    for (; plp; plp = plp->next) {
        pileup_t *p = plp;

        if (p->next)
            _mm_prefetch(p->next, _MM_HINT_T0);

        if (p->qual < opts->min_qual)
            continue;

        if (p->ref_skip)
            continue;

#ifdef K2
        int hb = 0;
#define _ 0
        static int X[16] = {_,0,1,_,2,_,_,_,3,_,_,_,_,_,_,_};
        int i, base1 = X[p->base4];
        for (i = p->seq_offset-K2; i <= p->seq_offset+K2; i++) {
            int base = i >= 0 && i < p->b.core.l_qseq ? X[bam_seqi(seq,i)] : _;
            hb = (hb<<2)|base;
        }
#undef _
#endif

        const bam1_t *b = &p->b;
        uint8_t base = p->base4;
        uint8_t *qual_arr = bam_get_qual(b);
        uint8_t qual = p->qual;
        //qual = qual*qual/40+1;
        if (qual == 255 || (qual == 0 && *qual_arr == 255))
            qual = default_qual;

#ifdef K2
        //qual = qual * hashN[hb] / hash1[base1];
        qual -= -TENOVERLOG10*log(hashN[hb] / (hash1[base1]+.1));
        if (qual < 1)
            qual = 1;
#endif

        // =ACM GRSV TWYH KDBN *
        static int L[32] = {
            5,0,1,5, 2,5,5,5, 3,5,5,5, 5,5,5,5,
            4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
        };

        // convert from sam base to acgt*n order.
        base = L[base];

        double MM, __, _M, oo, oM, o_, uu, um, mm, qe;

        // Correction for mapping quality.  Maybe speed up via lookups?
        // Cannot nullify mapping quality completely.  Lots of (true)
        // SNPs means low mapping quality.  (Ideally need to know
        // hamming distance to next best location.)

        if (flags & CONS_MQUAL) {
            int mqual = b->core.qual;
            if (opts->nm_adjust) {
                //mqual /= (nm_local(p, b, pos)+1);
                mqual /= (nm_local(p, b, b->core.pos + p->seq_offset+1)+1);
                mqual *= 1 + 2*(0.5-(td>30?30:td)/60.0); // depth fudge
            }

            // higher => call more; +FP, -FN
            // lower  => call less; -FP, +FN
            mqual *= opts->scale_mqual;

            // Drop these?  They don't seem to ever help.
            if (mqual < opts->low_mqual)
                mqual = opts->low_mqual;
            if (mqual > opts->high_mqual)
                mqual = opts->high_mqual;

            double _p = 1-q2p[qual];
            double _m = mqual_pow[mqual];
            qual = ph_log(1-(_m * _p + (1 - _m)/4)); // CURRENT
            //qual = ph_log(1-_p*_m); // testing
            //qual *= 6/sqrt(td);
        }

        /* Quality 0 should never be permitted as it breaks the maths */
        if (qual < 1)
            qual = 1;

        double poly = poly_len(p, b, b->core.pos + p->seq_offset+1);
#ifdef DO_POLY_DIST
        poly_dist[bam_is_rev(b)][MIN(99,(int)poly)]++;
#endif

        // EXPERIMENTAL
        // Adjust qual based on homopolymer length.
        // Affects different platforms by differing amounts.
        // May wish to further separate to qual2 and qual3 for ins and del?
        int qual2 = MAX(1, qual-(poly-2)*cp->poly_mul);

        /* MM=match  _M=half-match  __=mismatch */
        __ = cp->p__[qual];       // neither match
        MM = cp->pMM[qual] - __;  // both match
        _M = cp->p_M[qual] - __;  // one allele only (half match)

        /* observation ACGT, but against hypothesis ** or *base */
        oo = cp->poo[qual2] - __;
        oM = cp->poM[qual2] - __;
        o_ = cp->po_[qual2] - __;

        /* observation * */
        uu = cp->puu[qual2] - __;
        um = cp->pum[qual2] - __;
        mm = cp->pmm[qual2] - __;

        if (flags & CONS_DISCREP) {
            qe = q2p[qual];
            sumsC[base] += 1 - qe;
        }


        counts[base]++;
#ifdef DO_FRACT
        counts2[bam_is_rev(b)][base]++;
#endif

        // oM should never be higher than _M for actual bases!  or...
        //printf("base %d@%d MM %f _M %f oM %f\n", base, qual, MM, _M, oM);

        switch (base) {
        case 0: // A
            S[0]  += MM;
            S[1]  += _M;
            S[2]  += _M;
            S[3]  += _M;
            S[4]  += oM;
            S[8]  += o_;
            S[11] += o_;
            S[13] += o_;
            S[14] += oo;
            break;

        case 1: // C
            S[1]  += _M;
            S[5]  += MM;
            S[6]  += _M;
            S[7]  += _M;
            S[8]  += oM;
            S[4]  += o_;
            S[11] += o_;
            S[13] += o_;
            S[14] += oo;

            //fprintf(stderr, "%d %f %f %f\n", qual, MM+__, oo+__, MM-oo);
            break;

        case 2: // G
            S[ 2] += _M;
            S[ 6] += _M;
            S[ 9] += MM;
            S[10] += _M;
            S[11] += oM;
            S[4]  += o_;
            S[8]  += o_;
            S[13] += o_;
            S[14] += oo;
            break;

        case 3: // T
            S[ 3] += _M; // _m
            S[ 7] += _M;
            S[10] += _M;
            S[12] += MM; // mm
            S[13] += oM;
            S[4]  += o_;
            S[8]  += o_;
            S[11] += o_;
            S[14] += oo;
            // S[14] oo

            break;

        case 4: // *
            //   under       under       under       under   agree-no-base
            S[0] += uu; S[1 ]+= uu; S[2 ]+= uu; S[3 ]+= uu; S[4 ]+= um;
                        S[5 ]+= uu; S[6 ]+= uu; S[7 ]+= uu; S[8 ]+= um;
                                    S[9 ]+= uu; S[10]+= uu; S[11]+= um;
                                                S[12]+= uu; S[13]+= um;
                                                            S[14]+= mm;
            break;

        case 5: /* N => equal weight to all A,C,G,T but not a pad */
            S[0] += MM; S[1 ]+= MM; S[2 ]+= MM; S[3 ]+= MM; S[4 ]+= oM;
                        S[5 ]+= MM; S[6 ]+= MM; S[7 ]+= MM; S[8 ]+= oM;
                                    S[9 ]+= MM; S[10]+= MM; S[11]+= oM;
                                                S[12]+= MM; S[13]+= oM;
                                                            S[14]+= oo;
            break;
        }

        depth++;
    }

#ifdef DO_POLY_DIST
    // Or compute mean and s.d per strand.
    // Then compare likelihood of strands coming from the same distribution?
    // eg s.d=0.59 vs mean=3.41 sd=0.54... hmm
    //
    // Or compare ratio of most frequent to next most frequent, for each
    // strand.

    int d1 = 0, d2 = 0;
    double nd1 = 0, nd2 = 0;
    int k;
    for (k = 0; k < 100; k++) {
        if (!poly_dist[0][k] && !poly_dist[1][k])
            continue;

//        fprintf(stdout, "%ld %d %2d %2d\n", pos, k, poly_dist[0][k], poly_dist[1][k]);
        d1 += (k+1)*poly_dist[0][k];
        d2 += (k+1)*poly_dist[1][k];
        nd1 += poly_dist[0][k];
        nd2 += poly_dist[1][k];
    }
//    printf("Avg = %f / %f %f / %f / %f\n",
//           (d1+d2+1)/(nd1+nd2+1.),
//           (d1+1)/(nd1+1.), (d2+1)/(nd2+1.),
//           (d2+1)/(nd2+1.) - (d1+1)/(nd1+1.),
//           ((d2+1)/(nd2+1.) - (d1+1)/(nd1+1.)) / ((d1+d2+1)/(nd1+nd2+1.)));

    // Find the top two frequent lengths
    int n1 = 0, n2 = 0, l1 = 0, l2 = 0;
    for (k = 0; k < 100; k++) {
        int poly12 = poly_dist[0][k]+poly_dist[1][k];
        if (n1 < poly12) {
            n2 = n1; l2 = l1;
            n1 = poly12;
            l1 = k;
        } else if (n2 < poly12) {
            n2 = poly12;
            l2 = k;
        }
    }

    const double N = 5;
    nd1 += 1;
    nd2 += 1;

    // l1 is most common length
    int pn1p = poly_dist[0][l1];
    int pn1m = poly_dist[1][l1];
    // l2 2nd most common
    int pn2p = poly_dist[0][l2];
    int pn2m = poly_dist[1][l2];

    // ratio if two most common lengths on +
    double s1 = (pn1p+N) / (pn2p+N); s1 = s1>1?1/s1:s1;
    // ratio if two most common lengths on -
    double s2 = (pn1m+N) / (pn2m+N); s2 = s2>1?1/s2:s2;

    // ratio of s1 and s2 to identify strand bias
    double sbias = s1 / s2; sbias = sbias>1?1/sbias:sbias;

    if (pn2p+pn2m > 0 && l1 != l2) {
//        printf("len %d,%d  + %d,%d  - %d,%d\tbias = %f %f, %f %f\t%ld\n",
//               l1, l2,  pn1p, pn2p,   pn1m, pn2m,
//               s1, s2, sbias, sqrt(sbias)-1, pos);

        // adjust score for het indels
        // sbias is close to 0 for strong strand bias, and 1 for none
        sbias = 10*log(sbias);//+.5);
        S[ 4] += sbias; // A*
        S[ 8] += sbias; // C*
        S[11] += sbias; // G*
        S[13] += sbias; // T*
    } else {
        sbias = 0;
    }
#endif

    /* We've accumulated stats, so now we speculate on the consensus call */
    double shift, max, max_het, norm[15];
    int call = 0, het_call = 0;
    double tot1 = 0, tot2 = 0;

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

#ifdef DO_FRACT
    // Filter by --min-depth and --het-fract.
    // Also add a slight adjustment for strand bias.
    for (j = 0; j < 15; j++) {
        if (j == 0 || j == 5 || j == 9 || j == 12 || j == 14)
            continue;

        double c1p = counts2[0][map_het[j]%5];
        double c1m = counts2[1][map_het[j]%5];
        double c2p = counts2[0][map_het[j]/5];
        double c2m = counts2[1][map_het[j]/5];

        double c1 = c1p + c1m;
        double c2 = c2p + c2m;

        if (c1 && c2) {
            // Slight decrease in confidence if strong strand bias.
            const int N = 10; // avoid low sample size problems
            double b1 = 1 - (N+MIN(c1p,c1m))/(N+MAX(c1p,c1m));
            double b2 = 1 - (N+MIN(c2p,c2m))/(N+MAX(c2p,c2m));
            if (b1 > 0.5) S[j] -= b1;
            if (b2 > 0.5) S[j] -= b2;

            // Fraction based filtering, via --min-depth and --het-fract opts.
            c1 += 1e-5;
            c2 += 1e-5;
            if (c2 > c1) {
                double tmp = c2;
                c2 = c1;
                c1 = tmp;
            }

            if (c2 < opts->min_depth)
                S[j] -= 100;
            if (c2 / (c1+1e-5) <= opts->het_fract)
                S[j] -= 100;
        }
    }
#endif

#ifdef DO_HDW
    /*
     * Apply Hardy-Weinberg statistics for heterozygous sites.
     * This helps, but it also loses sensitivity a little.
     */
    for (j = 0; j < 15; j++) {
        if (j == 0 || j == 5 || j == 9 || j == 12 || j == 14)
            continue;

        double c1 = counts[map_het[j]%5];
        double c2 = counts[map_het[j]/5];

        if (c1 && c2) {
            c1 += 1e-5;
            c2 += 1e-5;
            if (c2 > c1) {
                double tmp = c2;
                c2 = c1;
                c1 = tmp;
            }

            // Limit depth for HW as we'll have an allele freq difference,
            // even if it's just caused by alignment reference bias.
            double c12 = c1+c2;
            if (c12 > 20) {
                c2 *= 20/(c12);
                c12 = 20;
                c1  = 20-c2;
            }

            // Helps a little, especially reducing FN deletions.
            c1+=1;
            c2+=1;
            c12+=2;
            S[j] += lnbinprobhalf(c12, c2) + fast_log2(c12)*0.69+.2;
        }
    }
#endif

    for (j = 0; j < 15; j++) {
        S[j] += cp->lprior15[j];
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

    for (j = 0; j < 15; j++) {
        norm[j]    += tot1;
        norm[14-j] += tot2;
        tot1 += S[j];
        tot2 += S[14-j];
    }

    /* And store result */
    if (!depth || depth == counts[5] /* all N */) {
        cons->call = 4; /* N */
        cons->het_call = 0;
        cons->het_logodd = 0;
        cons->phred = 0;
        cons->depth = 0;
        cons->discrep = 0;
        return 0;
    }

    cons->depth = depth;

    /* Call */
    if (norm[call] == 0) norm[call] = DBL_MIN;
    // Approximation of phred for when S[call] ~= 1 and norm[call]
    // is small.  Otherwise we need the full calculation.
    int ph;
    if (S[call] == 1 && norm[call] < .01)
        ph = ph_log(norm[call]) + .5;
    else
        ph = ph_log(1-S[call]/(norm[call]+S[call])) + .5;

    cons->call     = map_sing[call];
    cons->phred = ph < 0 ? 0 : ph;

    if (norm[het_call] == 0) norm[het_call] = DBL_MIN;
    ph = TENLOG2OVERLOG10 * (fast_log2(S[het_call])
                             - fast_log2(norm[het_call])) + .5;

    cons->het_call = map_het[het_call];
    cons->het_logodd = ph;

    /* Compute discrepancy score */
    if (flags & CONS_DISCREP) {
        double m = sumsC[0]+sumsC[1]+sumsC[2]+sumsC[3]+sumsC[4];
        double c;
        if (cons->het_logodd > 0)
            c = sumsC[cons->het_call%5] + sumsC[cons->het_call/5];
        else
            c = sumsC[cons->call];
        cons->discrep = (m-c)/sqrt(m);
    }

    return 0;
}

// If opts->gap5 is MODE_MIXED then we use two different parameter
// sets, favouring cp_p for precision and cp_r for recall.  Otherwise it's
// always cp_r only.
//
// When both calls equal, we return the same result.  When they differ,
// we adjust qual based on accurate vs recall profiles.
int calculate_consensus_gap5m(hts_pos_t pos, int flags, int depth,
                              pileup_t *plp, consensus_opts *opts,
                              consensus_t *cons, int default_qual,
                              cons_probs *cp_r, cons_probs *cp_p) {
    if (opts->mode != MODE_MIXED)
        return calculate_consensus_gap5(pos, flags, depth, plp, opts,
                                        cons, default_qual,
                                        opts->mode == MODE_PRECISE
                                            ? cp_p : cp_r);

    // EXPERIMENTAL: mixed mode
    consensus_t consP, consR;
    // Favours precision
    calculate_consensus_gap5(pos, flags, depth, plp, opts,
                             &consP, default_qual, cp_p);
    // Favours recall
    calculate_consensus_gap5(pos, flags, depth, plp, opts,
                             &consR, default_qual, cp_r);

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

    // Initial starting point is precise mode
    memcpy(cons, &consP, sizeof(consP));

    if (consP.phred > 0 && consR.phred > 0 && consP.call == consR.call) {
        // Both strategies match as HOM
        // Boost qual as both in agreement
        cons->phred += MIN(20, consR.phred);

    } else if (consP.het_logodd >= 0 && consR.het_logodd >= 0 &&
               consP.het_call == consR.het_call) {
        // Both strategies match as HET
        // Boost qual as both in agreement
        cons->het_logodd += MIN(20, consR.het_logodd);

    } else if (consP.het_logodd >= 0) {
        // Accurate method claims heterozygous, so go with it.
        // However sensitive method disagrees, so reduce qual a little.
        int q2 = MAX(consR.phred, consR.het_logodd);
        cons->het_logodd = MAX(1, (cons->het_logodd - q2/2));

    } else if (consR.het_logodd >= 70) {
        // Accurate is homozygous and consR is het, so we go with it instead
        // but at a lower quality value.
        // TODO: may wish to check HET is consistent with HOM? Very unlikely
        // not to be though.
        int q1 = consP.phred;
        int q2 = consR.het_logodd;
        memcpy(cons, &consR, sizeof(consR));
        cons->het_logodd = MIN(15, MAX((q2-q1*2)/2, 1+q2/(q1+1.0)));

    } else if (consR.het_logodd >= 0) {
        // As above, but low quality
        int q1 = consP.phred;
        int q2 = consR.het_logodd;
        memcpy(cons, &consR, sizeof(consR));
        cons->het_logodd = MAX(1,q2 - 0.3*q1)
            + 5*(consP.het_call == consR.het_call);
        cons->phred = 0;

    } else if (consR.het_logodd < 0) {
        // Neither are heterozygous, but differing in phred call (V rare).
        // Pick highest qual, after some scaling?
        consR.phred = consR.phred / 2;
        if (consR.phred > consP.phred)
            memcpy(cons, &consR, sizeof(consR));
        cons->phred = MAX(10, cons->phred);
    }

    return 0;
}

/* --------------------------------------------------------------------------
 * Main processing logic
 */

static void dump_fastq(consensus_opts *opts,
                       const char *name,
                       const char *seq, size_t seq_l,
                       const char *qual, size_t qual_l) {
    enum format fmt = opts->fmt;
    int line_len = opts->line_len;
    FILE *fp = opts->fp_out;

    fprintf(fp, "%c%s\n", ">@"[fmt==FASTQ], name);
    size_t i;
    for (i = 0; i < seq_l; i += line_len)
        fprintf(fp, "%.*s\n", (int)MIN(line_len, seq_l - i), seq+i);

    if (fmt == FASTQ) {
        fprintf(fp, "+\n");
        for (i = 0; i < seq_l; i += line_len)
            fprintf(fp, "%.*s\n", (int)MIN(line_len, seq_l - i), qual+i);
    }
}

//---------------------------------------------------------------------------

/*
 * Reads a single alignment record, using either the iterator
 * or a direct sam_read1 call.
 */
static int readaln2(void *dat, samFile *fp, sam_hdr_t *h, bam1_t *b) {
    consensus_opts *opts = (consensus_opts *)dat;

    for (;;) {
        int ret = opts->iter
            ? sam_itr_next(fp, opts->iter, b)
            : sam_read1(fp, h, b);
        if (ret < 0)
            return ret;

        // Apply hard filters
        if (opts->incl_flags && !(b->core.flag & opts->incl_flags))
            continue;
        if (opts->excl_flags &&  (b->core.flag & opts->excl_flags))
            continue;
        if (b->core.qual < opts->min_mqual)
            continue;

        return ret;
    }
}

/* --------------------------------------------------------------------------
 * A simple summing algorithm, either pure base frequency, or by
 * weighting them according to their quality values.
 *
 * This is crude, but easy to understand and fits with several
 * standard pileup criteria (eg COG-UK / CLIMB Covid-19 seq project).
 *
 *
 * call1 / score1 is the highest scoring allele.
 * call2 / score2 is the second highest scoring allele.
 *
 * Het_fract:  score2/score1
 * Call_fract: score1 or score1+score2 over total score
 * Min_depth:  minimum total depth of unfiltered bases (above qual/mqual)
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
static int calculate_consensus_simple(const pileup_t *plp,
                                      consensus_opts *opts, int *qual) {
    int i, min_qual = opts->min_qual;
    int tot_depth = 0;

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

    // Accumulate
    for (; plp; plp = plp->next) {
        const pileup_t *p = plp;
        if (p->next)
            _mm_prefetch(p->next, _MM_HINT_T0);

        int q = p->qual;
        if (q < min_qual)
            // Should we still record these in freq[] somewhere so
            // we can use them in the fracts?
            // Difference between >= X% of high-qual bases calling Y
            // and >= X% of all bases are high-quality Y calls.
            continue;

        //int b = p->is_del ? 16 : bam_seqi(bam_get_seq(&p->b), p->seq_offset);
        int b = p->base4;

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
        tot_depth++;
    }

    // Total usable depth
    int tscore = 0;
    for (i = 0; i < 5; i++)
        tscore += score[1<<i];

    // Best and second best potential calls
    int call1  = 15, call2 = 15;
    int score1 = 0,  score2 = 0;
    for (i = 0; i < 5; i++) {
        int c = 1<<i; // A C G T *
        if (score1 < score[c]) {
            score2 = score1;
            call2  = call1;
            score1 = score[c];
            call1  = c;
        } else if (score2 < score[c]) {
            score2 = score[c];
            call2  = c;
        }
    }

    // Work out which best and second best are usable as a call
    int used_score = score1;
    int used_base  = call1;
    if (score2 >= opts->het_fract * score1 && opts->ambig) {
        used_base  |= call2;
        used_score += score2;
    }

    // N is too shallow, or insufficient proportion of total
    if (tot_depth  < opts->min_depth ||
        used_score < opts->call_fract * tscore) {
        // But note shallow gaps are still called gaps, not N, as
        // we're still more confident there is no base than it is
        // A, C, G or T.
        used_base = call1 == 16 ? 16 : 0; // * or N
    }

    // Our final call.  "?" shouldn't be possible to generate
    const char *het =
        "NACMGRSVTWYHKDBN"
        "*ac?g???t???????";

    //printf("%c %d\n", het[used_base], tot_depth);
    if (qual)
        *qual = used_base ? 100.0 * used_score / tscore : 0;

    return het[used_base];
}

static int empty_pileup2(consensus_opts *opts, sam_hdr_t *h, int tid,
                         hts_pos_t start, hts_pos_t end) {
    const char *name = sam_hdr_tid2name(h, tid);
    hts_pos_t i;

    int err = 0;
    for (i = start; i < end; i++)
        err |= fprintf(opts->fp_out, "%s\t%"PRIhts_pos"\t0\t0\tN\t0\t*\t*\n", name, i+1) < 0;

    return err ? -1 : 0;
}

/*
 * Returns 0 on success
 *        -1 on failure
 */
static int basic_pileup(void *cd, samFile *fp, sam_hdr_t *h, pileup_t *p,
                        int depth, hts_pos_t pos, int nth, int is_insert) {
    unsigned char *qp, *cp;
    char *rp;
    int ref, cb, cq;
    consensus_opts *opts = (consensus_opts *)cd;
    int tid = p->b.core.tid;

//    opts->show_ins=0;
//    opts->show_del=1;
    if (!opts->show_ins && nth)
        return 0;

    if (opts->iter) {
        if (opts->iter->beg >= pos || opts->iter->end < pos)
            return 0;
    }

    if (opts->all_bases) {
        if (tid != opts->last_tid && opts->last_tid >= 0) {
            hts_pos_t len = sam_hdr_tid2len(opts->h, opts->last_tid);
            if (opts->iter)
                len =  MIN(opts->iter->end, len);
            if (empty_pileup2(opts, opts->h, opts->last_tid, opts->last_pos,
                              len) < 0)
                return -1;
            if (tid >= 0) {
                if (empty_pileup2(opts, opts->h, tid,
                                  opts->iter ? opts->iter->beg : 0,
                                  pos-1) < 0)
                    return -1;
            }
        }
        if (opts->last_pos >= 0 && pos > opts->last_pos+1) {
            if (empty_pileup2(opts, opts->h, p->b.core.tid, opts->last_pos,
                              pos-1) < 0)
                return -1;
        } else if (opts->last_pos < 0) {
            if (empty_pileup2(opts, opts->h, p->b.core.tid,
                              opts->iter ? opts->iter->beg : 0, pos-1) < 0)
                return -1;
        }
    }

    if (opts->mode != MODE_SIMPLE) {
        consensus_t cons;
        calculate_consensus_gap5m(pos, opts->use_mqual ? CONS_MQUAL : 0,
                                  depth, p, opts, &cons, opts->default_qual,
                                  &cons_prob_recall, &cons_prob_precise);
        if (cons.het_logodd > 0 && opts->ambig) {
            cb = "AMRWa" // 5x5 matrix with ACGT* per row / col
                 "MCSYc"
                 "RSGKg"
                 "WYKTt"
                 "acgt*"[cons.het_call];
            cq = cons.het_logodd;
        } else{
            cb = "ACGT*"[cons.call];
            cq = cons.phred;
        }
        if (cq < opts->cons_cutoff && cb != '*') {
            cb = 'N';
            cq = 0;
        }
    } else {
        cb = calculate_consensus_simple(p, opts, &cq);
    }
    if (cb < 0)
        return -1;

    if (!p)
        return 0;

    if (!opts->show_del && cb == '*')
        return 0;

    /* Ref, pos, nth, score, seq, qual */
    kstring_t *ks = &opts->ks_line;
    ks->l = 0;
    ref = p->b.core.tid;
    rp = (char *)sam_hdr_tid2name(h, ref);

    int err = 0;
    err |= kputs(rp, ks)    < 0;
    err |= kputc_('\t', ks) < 0;
    err |= kputw(pos, ks)   < 0;
    err |= kputc_('\t', ks) < 0;
    err |= kputw(nth, ks)   < 0;
    err |= kputc_('\t', ks) < 0;
    err |= kputw(depth, ks) < 0;
    err |= kputc_('\t', ks) < 0;
    err |= kputc_(cb, ks)   < 0;
    err |= kputc_('\t', ks) < 0;
    err |= kputw(cq, ks)    < 0;
    err |= kputc_('\t', ks) < 0;
    if (err)
        return -1;

    /* Seq + qual at predetermined offsets */
    if (ks_resize(ks, ks->l + depth*2 + 2) < 0)
        return -1;

    cp = (unsigned char *)ks->s + ks->l;
    ks->l += depth*2 + 2;
    qp = cp+depth+1;
    for (; p; p = p->next) {
        // Too tight a loop to help much, but some benefit still
        if (p->next && p->next->next)
            _mm_prefetch(p->next->next, _MM_HINT_T0);
        if (p->b_is_rev) {
            *cp++ = p->base == '*' ? '#' : tolower(p->base);
        } else {
            *cp++ = p->base;
        }
        *qp++ = MIN(p->qual,93) + '!';
    }
    *cp++ = '\t';
    *qp++ = '\n';
    if (fwrite(ks->s, 1, ks->l, opts->fp_out) != ks->l)
        return -1;

    opts->last_pos = pos;
    opts->last_tid = tid;

    return 0;
}

static int basic_fasta(void *cd, samFile *fp, sam_hdr_t *h, pileup_t *p,
                       int depth, hts_pos_t pos, int nth, int is_insert) {
    int cb, cq;
    consensus_opts *opts = (consensus_opts *)cd;
    int tid = p->b.core.tid;
    kstring_t *seq  = &opts->ks_ins_seq;
    kstring_t *qual = &opts->ks_ins_qual;

    if (!opts->show_ins && nth)
        return 0;

    if (opts->iter) {
        if (opts->iter->beg >= pos || opts->iter->end < pos)
            return 0;
    }

    if (tid != opts->last_tid) {
        if (opts->last_tid != -1) {
            if (opts->all_bases) {
                int i, N;
                if (opts->iter) {
                    opts->last_pos = MAX(opts->last_pos, opts->iter->beg-1);
                    N = opts->iter->end;
                } else {
                    N = INT_MAX;
                }
                N = MIN(N, sam_hdr_tid2len(opts->h, opts->last_tid))
                    - opts->last_pos;
                if (N > 0) {
                    if (ks_expand(seq, N+1) < 0)
                        return -1;
                    if (ks_expand(qual, N+1) < 0)
                        return -1;
                    for (i = 0; i < N; i++) {
                        seq->s[seq->l++] = 'N';
                        qual->s[qual->l++] = '!';
                    }
                    seq->s[seq->l] = 0;
                    qual->s[qual->l] = 0;
                }
            }
            dump_fastq(opts, sam_hdr_tid2name(opts->h, opts->last_tid),
                       seq->s, seq->l, qual->s, qual->l);
        }

        seq->l = 0; qual->l = 0;
        opts->last_tid = tid;
//        if (opts->all_bases)
//            opts->last_pos = 0;
        if (opts->iter)
            opts->last_pos = opts->iter->beg;
        else
            opts->last_pos = opts->all_bases ? 0 : pos-1;
    }

    // share this with basic_pileup
    if (opts->mode != MODE_SIMPLE) {
        consensus_t cons;
        calculate_consensus_gap5m(pos, opts->use_mqual ? CONS_MQUAL : 0,
                                  depth, p, opts, &cons, opts->default_qual,
                                  &cons_prob_recall, &cons_prob_precise);
        if (cons.het_logodd > 0 && opts->ambig) {
            cb = "AMRWa" // 5x5 matrix with ACGT* per row / col
                 "MCSYc"
                 "RSGKg"
                 "WYKTt"
                 "acgt*"[cons.het_call];
            cq = cons.het_logodd;
        } else{
            cb = "ACGT*"[cons.call];
            cq = cons.phred;
        }
        if (cq < opts->cons_cutoff && cb != '*' &&
            cons.het_call % 5 != 4 && cons.het_call / 5 != 4) {
            // het base/* keeps base or * as most likely pure call, else N.
            // This is because we don't have a traditional way of representing
            // base or not-base ambiguity.
            cb = 'N';
            cq = 0;
        }
    } else {
        cb = calculate_consensus_simple(p, opts, &cq);
    }
    if (cb < 0)
        return -1;

    if (!p)
        return 0;

    if (!opts->show_del && cb == '*') {
        opts->last_pos = pos;
        opts->last_tid = tid;
        return 0;
    }
    if (opts->mark_ins && nth && cb != '*') {
        kputc('_', seq);
        kputc('_', qual);
    }

    // end of share

    // Append consensus base/qual to seqs
    if (pos > opts->last_pos) {
        if (opts->last_pos >= 0 || opts->all_bases) {
            // FIXME: don't expand qual if fasta
            if (ks_expand(seq,  pos - opts->last_pos) < 0 ||
                ks_expand(qual, pos - opts->last_pos) < 0)
                return -1;
            memset(seq->s  + seq->l,  'N', pos - (opts->last_pos+1));
            memset(qual->s + qual->l, '!', pos - (opts->last_pos+1));
            seq->l  += pos - (opts->last_pos+1);
            qual->l += pos - (opts->last_pos+1);
        }
    }
    if ((nth && opts->show_ins && cb != '*')
        || cb != '*' || (pos > opts->last_pos && opts->show_del)) {
        int err = 0;
        err |= kputc(cb, seq) < 0;
        err |= kputc(MIN(cq, '~'-'!')+'!', qual) < 0;
        if (err)
            return -1;
    }

    opts->last_pos = pos;
    opts->last_tid = tid;

    return 0;
}

// END OF NEW PILEUP
//---------------------------------------------------------------------------

static void usage_exit(FILE *fp, int exit_status) {
    fprintf(fp, "Usage: samtools consensus [options] <in.bam>\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "  -r, --region REG      Limit query to REG. Requires an index\n");
    fprintf(fp, "  -f, --format FMT      Output in format FASTA, FASTQ or PILEUP [FASTA]\n");
    fprintf(fp, "  -l, --line-len INT    Wrap FASTA/Q at line length INT [70]\n");
    fprintf(fp, "  -o, --output FILE     Output consensus to FILE\n");
    fprintf(fp, "  -m, --mode STR        Switch consensus mode to \"simple\"/\"bayesian\" [bayesian]\n");
    fprintf(fp, "  -a                    Output all bases (start/end of reference)\n");
    fprintf(fp, "  --rf, --incl-flags STR|INT\n");
    fprintf(fp, "                        Only include reads with any flag bit set [0]\n");
    fprintf(fp, "  --ff, --excl-flags STR|INT\n");
    fprintf(fp, "                        Exclude reads with any flag bit set\n");
    fprintf(fp, "                        [UNMAP,SECONDARY,QCFAIL,DUP]\n");
    fprintf(fp, "  --min-MQ INT          Exclude reads with mapping quality below INT [0]\n");
    fprintf(fp, "  --min-BQ INT          Exclude reads with base quality below INT [0]\n");
    fprintf(fp, "  --show-del yes/no     Whether to show deletion as \"*\" [no]\n");
    fprintf(fp, "  --show-ins yes/no     Whether to show insertions [yes]\n");
    fprintf(fp, "  --mark-ins            Add '+' before every inserted base/qual [off]\n");
    fprintf(fp, "  -A, --ambig           Enable IUPAC ambiguity codes [off]\n");
    fprintf(fp, "\nFor simple consensus mode:\n");
    fprintf(fp, "  -q, --(no-)use-qual   Use quality values in calculation [off]\n");
    fprintf(fp, "  -c, --call-fract INT  At least INT portion of bases must agree [0.75]\n");
    fprintf(fp, "  -d, --min-depth INT   Minimum depth of INT [2]\n");
    fprintf(fp, "  -H, --het-fract INT   Minimum fraction of 2nd-most to most common base [0.15]\n");
    fprintf(fp, "\nFor default \"Bayesian\" consensus mode:\n");
    fprintf(fp, "  -C, --cutoff C        Consensus cutoff quality C [10]\n");
    fprintf(fp, "      --(no-)adj-qual   Modify quality with local minima [on]\n");
    fprintf(fp, "      --(no-)use-MQ     Use mapping quality in calculation [on]\n");
    fprintf(fp, "      --(no-)adj-MQ     Modify mapping quality by local NM [on]\n");
    fprintf(fp, "      --NM-halo INT     Size of window for NM count in --adj-MQ [50]\n");
    fprintf(fp, "      --scale-MQ FLOAT  Scale mapping quality by FLOAT [1.00]\n");
    fprintf(fp, "      --low-MQ  INT     Cap minimum mapping quality [1]\n");
    fprintf(fp, "      --high-MQ INT     Cap maximum mapping quality [60]\n");
    fprintf(fp, "      --P-het FLOAT     Probability of heterozygous site[%.1e]\n",
            P_HET);
    fprintf(fp, "      --P-indel FLOAT   Probability of indel sites[%.1e]\n",
            P_INDEL);
    fprintf(fp, "      --het-scale FLOAT Heterozygous SNP probability multiplier[%.1e]\n",
            P_HET_SCALE);
    fprintf(fp, "  -p, --homopoly-fix    Spread low-qual bases to both ends of homopolymers\n");
    fprintf(fp, "      --homopoly-score FLOAT\n"
                "                        Qual fraction adjustment for -p option [%g]\n", P_HOMOPOLY);
    fprintf(fp, "  -t, --qual-calibration FILE / :config (see man page)\n");
    fprintf(fp, "                        Load quality calibration file\n");
    fprintf(fp, "\n");
    fprintf(fp, "  -X, --config STR      Use pre-defined configuration set. STR from:\n");
    fprintf(fp, "                        hiseq, hifi, r10.4_sup, r10.4_dup and ultima\n");

    fprintf(fp, "\nGlobal options:\n");
    sam_global_opt_help(fp, "-.---@-.");
    exit(exit_status);
}

int main_consensus(int argc, char **argv) {
    int c, ret = 1;

    consensus_opts opts = {
        // User options
        .mode         = MODE_RECALL,
        .use_qual     = 0,
        .min_qual     = 0,
        .adj_qual     = 1,
        .use_mqual    = 1,
        .scale_mqual  = 1.00,
        .nm_adjust    = 1,
        .nm_halo      = 50,
        .sc_cost      = 60,
        .low_mqual    = 1,
        .high_mqual   = 60,
        .min_depth    = 1,
        .call_fract   = 0.75,
        .het_fract    = 0.5,
        .het_only     = 0,
        .fmt          = FASTA,
        .cons_cutoff  = 10,
        .ambig        = 0,
        .line_len     = 70,
        .default_qual = 10,
        .all_bases    = 0,
        .show_del     = 0,
        .show_ins     = 1,
        .mark_ins     = 0,
        .incl_flags   = 0,
        .excl_flags   = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP,
        .min_mqual    = 0,
        .P_het        = P_HET,
        .P_indel      = P_INDEL,
        .het_scale    = P_HET_SCALE,
        .homopoly_fix = 0,
        .homopoly_redux = 0.01,

        // Internal state
        .ks_line      = {0,0},
        .ks_ins_seq   = {0,0},
        .ks_ins_qual  = {0,0},
        .fp           = NULL,
        .fp_out       = stdout,
        .iter         = NULL,
        .idx          = NULL,
        .last_tid     = -1,
        .last_pos     = -1,
    };

    set_qcal(&opts.qcal, QCAL_FLAT);

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', '-', '-', '@'),
        {"use-qual",           no_argument,       NULL, 'q'},
        {"no-use-qual",        no_argument,       NULL, 'q'+1000},
        {"adj-qual",           no_argument,       NULL, 'q'+100},
        {"no-adj-qual",        no_argument,       NULL, 'q'+101},
        {"use-MQ",             no_argument,       NULL, 'm'+1000},
        {"no-use-MQ",          no_argument,       NULL, 'm'+1001},
        {"adj-MQ",             no_argument,       NULL, 'm'+100},
        {"no-adj-MQ",          no_argument,       NULL, 'm'+101},
        {"NM-halo",            required_argument, NULL, 'h'+100},
        {"SC-cost",            required_argument, NULL, 'h'+101},
        {"scale-MQ",           required_argument, NULL, 14},
        {"low-MQ"   ,          required_argument, NULL,  9},
        {"high-MQ",            required_argument, NULL, 10},
        {"min-depth",          required_argument, NULL, 'd'},
        {"call-fract",         required_argument, NULL, 'c'},
        {"het-fract",          required_argument, NULL, 'H'},
        {"region",             required_argument, NULL, 'r'},
        {"format",             required_argument, NULL, 'f'},
        {"cutoff",             required_argument, NULL, 'C'},
        {"ambig",              no_argument,       NULL, 'A'},
        {"line-len",           required_argument, NULL, 'l'},
        {"default-qual",       required_argument, NULL, 1},
        {"het-only",           no_argument,       NULL, 6},
        {"show-del",           required_argument, NULL, 7},
        {"show-ins",           required_argument, NULL, 8},
        {"mark-ins",           no_argument,       NULL, 18},
        {"output",             required_argument, NULL, 'o'},
        {"incl-flags",         required_argument, NULL, 11},
        {"rf",                 required_argument, NULL, 11},
        {"excl-flags",         required_argument, NULL, 12},
        {"ff",                 required_argument, NULL, 12},
        {"min-MQ",             required_argument, NULL, 13},
        {"min-BQ",             required_argument, NULL, 16},
        {"P-het",              required_argument, NULL, 15},
        {"P-indel",            required_argument, NULL, 17},
        {"het-scale",          required_argument, NULL, 19},
        {"mode",               required_argument, NULL, 'm'},
        {"homopoly-fix",       no_argument,       NULL, 'p'},
        {"homopoly-score",     required_argument, NULL, 'p'+100},
        {"homopoly-redux",     required_argument, NULL, 'p'+200},
        {"qual-calibration",   required_argument, NULL, 't'},
        {"config",             required_argument, NULL, 'X'},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "@:qd:c:H:r:5f:C:aAl:o:m:pt:X:",
                            lopts, NULL)) >= 0) {
        switch (c) {
        case 'a': opts.all_bases++; break;
        case 'q': opts.use_qual=1; break;
        case 'q'+1000: opts.use_qual=0; break;
        case 'm'+1000: opts.use_mqual=1; break;
        case 'm'+1001: opts.use_mqual=0; break;
        case 14:  opts.scale_mqual = atof(optarg); break;
        case  9:  opts.low_mqual = atoi(optarg); break;
        case 10:  opts.high_mqual = atoi(optarg); break;
        case 'd': opts.min_depth = atoi(optarg); break;
        case 'c': opts.call_fract = atof(optarg); break;
        case 'H': opts.het_fract = atof(optarg); break;
        case 'r': opts.reg = optarg; break;
        case 'C': opts.cons_cutoff = atoi(optarg); break;
        case 'A': opts.ambig = 1; break;
        case 'p': opts.homopoly_fix = P_HOMOPOLY; break;
        case 'p'+100: opts.homopoly_fix = atof(optarg); break;
        case 'p'+200:
          // EXPERIMENTAL
          opts.homopoly_redux = atof(optarg); break;
        case 1:   opts.default_qual = atoi(optarg); break;
        case 6:   opts.het_only = 1; break;
        case 7:   opts.show_del = (*optarg == 'y' || *optarg == 'Y'); break;
        case 8:   opts.show_ins = (*optarg == 'y' || *optarg == 'Y'); break;
        case 18:  opts.mark_ins = 1; break;
        case 13:  opts.min_mqual = atoi(optarg); break;
        case 16:  opts.min_qual  = atoi(optarg); break;
        case 15:  opts.P_het = atof(optarg); break;
        case 17:  opts.P_indel = atof(optarg); break;
        case 19:  opts.het_scale = atof(optarg); break;
        case 'q'+100: opts.adj_qual = 1; break;
        case 'q'+101: opts.adj_qual = 0; break;
        case 'm'+100: opts.nm_adjust = 1; break;
        case 'm'+101: opts.nm_adjust = 0; break;
        case 'h'+100: opts.nm_halo = atoi(optarg); break;
        case 'h'+101: opts.sc_cost = atoi(optarg); break;

        case 'm': // mode
            if (strcasecmp(optarg, "simple") == 0) {
                opts.mode = MODE_SIMPLE;
            } else if (strcasecmp(optarg, "bayesian_m") == 0) {
                // EXPERIMENTAL:
                // A mixture of modified precise/recall params and a
                // blending of the two.  Sometimes helps a bit.
                opts.mode = MODE_MIXED;
            } else if (strcasecmp(optarg, "bayesian_p") == 0) {
                // EXPERIMENTAL:
                // favours precision
                opts.mode = MODE_PRECISE;
            } else if (strcasecmp(optarg, "bayesian_r") == 0 ||
                       strcasecmp(optarg, "bayesian") == 0) {
                // favours recall; the default
                opts.mode = MODE_RECALL;
            } else if (strcasecmp(optarg, "bayesian_116") == 0) {
                opts.mode = MODE_BAYES_116;
            } else {
                fprintf(stderr, "Unknown mode %s\n", optarg);
                return 1;
            }
            break;

        case 'l':
            if ((opts.line_len = atoi(optarg)) <= 0)
                opts.line_len = INT_MAX;
            break;

        case 'f':
            if (strcasecmp(optarg, "fasta") == 0) {
                opts.fmt = FASTA;
            } else if (strcasecmp(optarg, "fastq") == 0) {
                opts.fmt = FASTQ;
            } else if (strcasecmp(optarg, "pileup") == 0) {
                opts.fmt = PILEUP;
            } else {
                fprintf(stderr, "Unknown format %s\n", optarg);
                return 1;
            }
            break;

        case 'o':
            if (!(opts.fp_out = fopen(optarg, "w"))) {
                perror(optarg);
                return 1;
            }
            break;

        case 'X':
            if (strcasecmp(optarg, "hifi") == 0) {
                set_qcal(&opts.qcal, QCAL_HIFI);
                opts.mode = MODE_RECALL;
                opts.homopoly_fix = 0.3;
                opts.homopoly_redux = 0.01;
                opts.low_mqual = 5;
                opts.scale_mqual = 1.5;
                opts.het_scale = 0.37;
            } else if (strcasecmp(optarg, "hiseq") == 0) {
                opts.mode = MODE_RECALL;
                set_qcal(&opts.qcal, QCAL_HISEQ);
                opts.homopoly_redux = 0.01;
            } else if (strcasecmp(optarg, "r10.4_sup") == 0) {
                // Same as HiFi params, but ONT calibration table.
                // At higher depth, hifi params work well for ONT
                // when combined with ONT calibration chart.
                //
                // At lower depth we gain a bit from increasing homopoly_redux
                set_qcal(&opts.qcal, QCAL_ONT_R10_4_SUP);
                opts.mode = MODE_RECALL;
                opts.homopoly_fix = 0.3;
                opts.homopoly_redux = 0.01;
                opts.low_mqual = 5;
                opts.scale_mqual = 1.5;
                opts.het_scale = 0.37;

                // Also consider, for lower depth:
                // opts.homopoly_redux = 1;
                // opts.scale_mqual = 1;
                // opts.het_scale = 0.45;
            } else if (strcasecmp(optarg, "r10.4_dup") == 0) {
                // Just a copy of of HiFi for duplex currently until
                // we get a good truth set for calibration.
                set_qcal(&opts.qcal, QCAL_ONT_R10_4_DUP);
                opts.mode = MODE_RECALL;
                opts.homopoly_fix = 0.3;
                opts.homopoly_redux = 0.01;
                opts.low_mqual = 5;
                opts.scale_mqual = 1.5;
                opts.het_scale = 0.37;
            } else if (strcasecmp(optarg, "ultima") == 0) {
                // Very similar to HiFi, but with own calibration table
                opts.mode = MODE_RECALL;
                set_qcal(&opts.qcal, QCAL_ULTIMA);
                opts.homopoly_fix = 0.3;
                opts.homopoly_redux = 0.01;
                opts.het_scale = 0.37;
                opts.scale_mqual = 2;
                opts.low_mqual = 10;
            } else {
                // NB consider defaults that are a mixture of all above.
                // Options are all similar for all bar Illumina.
                // Unsure what :flat calibration table does to each of
                // these though.
                fprintf(stderr, "Unrecognised configuration name: \"%s\"\n",
                        optarg);
                return 1;
            }
            break;

        case 11:
            if ((opts.incl_flags = bam_str2flag(optarg)) < 0) {
                print_error("consensus", "could not parse --rf %s", optarg);
                return 1;
            }
            break;
        case 12:
            if ((opts.excl_flags = bam_str2flag(optarg)) < 0) {
                print_error("consensus", "could not parse --ff %s", optarg);
                return 1;
            }
            break;

        case 't': // --qual-calibration
            if (load_qcal(&opts.qcal, optarg) < 0) {
                print_error("consensus",
                            "failed to load quality calibration '%s'",
                            optarg);
                return -1;
            }
            break;

         default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

#if 0
    // Dump out the qcal table.  Useful for copying into the code above.
    int i;
    qcal_t *q = &opts.qcal;
    fprintf(stderr, "{");
    for (i = 0; i < 100; i++)
        fprintf(stderr, "%2d,%s", q->smap[i],(i+1)%10?" ":"\n");
    fprintf(stderr, "},\n{");
    for (i = 0; i < 100; i++)
        fprintf(stderr, "%2d,%s", q->umap[i],(i+1)%10?" ":"\n");
    fprintf(stderr, "},\n{");
    for (i = 0; i < 100; i++)
        fprintf(stderr, "%2d,%s", q->omap[i],(i+1)%10?" ":"\n");
    fprintf(stderr, "}\n");
#endif

    if (opts.mode != MODE_SIMPLE) {
        if (opts.mode == MODE_PRECISE)
            // More accuracy / precision, but a significant drop
            // in recall.
            consensus_init(opts.P_het, opts.P_indel,
                           0.3 * opts.het_scale, opts.homopoly_redux,
                           &opts.qcal, MODE_PRECISE, &cons_prob_precise);

        if (opts.mode == MODE_MIXED)
            // Blend these in when running in mixed mode, so we can
            // keep sensitivity but have a better joint quality to
            // reduce the FP rate.
            consensus_init(pow(opts.P_het, 0.7), pow(opts.P_indel, 0.7),
                           0.3 * opts.het_scale, opts.homopoly_redux,
                           &opts.qcal, MODE_PRECISE, &cons_prob_precise);

        // Better recall, at a cost of some accuracy (false positives)
        consensus_init(opts.P_het, opts.P_indel, opts.het_scale,
                       opts.mode == MODE_RECALL ? opts.homopoly_redux : 0.01,
                       &opts.qcal, MODE_RECALL, &cons_prob_recall);
    }

    if (argc != optind+1) {
        if (argc == optind) usage_exit(stdout, EXIT_SUCCESS);
        else usage_exit(stderr, EXIT_FAILURE);
    }
    opts.fp = sam_open_format(argv[optind], "r", &ga.in);
    if (opts.fp == NULL) {
        print_error_errno("consensus", "Cannot open input file \"%s\"",
                          argv[optind]);
        goto err;
    }
    if (ga.nthreads > 0)
        hts_set_threads(opts.fp, ga.nthreads);

    if (hts_set_opt(opts.fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        goto err;
    }

    if (!(opts.h = sam_hdr_read(opts.fp))) {
        fprintf(stderr, "Failed to read header for \"%s\"\n", argv[optind]);
        goto err;
    }

    if (opts.reg) {
        opts.idx = sam_index_load(opts.fp, argv[optind]);
        if (!opts.idx) {
            print_error("consensus", "Cannot load index for input file \"%s\"",
                        argv[optind]);
            goto err;
        }
        opts.iter = sam_itr_querys(opts.idx, opts.h, opts.reg);
        if (!opts.iter) {
            print_error("consensus", "Failed to parse region \"%s\"",
                        opts.reg);
            goto err;
        }
    }

    if (opts.fmt == PILEUP) {
        if (pileup_loop(opts.fp, opts.h, readaln2,
                        opts.mode != MODE_SIMPLE ? nm_init : NULL,
                        basic_pileup,
                        opts.mode != MODE_SIMPLE ? nm_free : NULL,
                        &opts) < 0)
            goto err;

        if (opts.all_bases) {
            int tid = opts.iter ? opts.iter->tid : opts.last_tid;
            int len = sam_hdr_tid2len(opts.h, tid);
            int pos = opts.last_pos;
            if (opts.iter) {
                len = MIN(opts.iter->end, len);
                pos = MAX(opts.iter->beg, pos);
            }
            if (empty_pileup2(&opts, opts.h, tid, pos, len) < 0)
                goto err;
        }
    } else {
        if (pileup_loop(opts.fp, opts.h, readaln2,
                        opts.mode != MODE_SIMPLE ? nm_init : NULL,
                        basic_fasta,
                        opts.mode != MODE_SIMPLE ? nm_free : NULL,
                        &opts) < 0)
            goto err;
        if (opts.all_bases) {
            // fill out terminator
            int tid = opts.iter ? opts.iter->tid : opts.last_tid;
            int len = sam_hdr_tid2len(opts.h, tid);
            int pos = opts.last_pos;
            if (opts.iter) {
                len = MIN(opts.iter->end, len);
                pos = MAX(opts.iter->beg, pos);
                opts.last_tid = opts.iter->tid;
            }
            if (pos < len) {
                if (ks_expand(&opts.ks_ins_seq,  len-pos+1) < 0)
                    goto err;
                if (ks_expand(&opts.ks_ins_qual, len-pos+1) < 0)
                    goto err;
                while (pos++ < len) {
                    opts.ks_ins_seq.s [opts.ks_ins_seq.l++] = 'N';
                    opts.ks_ins_qual.s[opts.ks_ins_qual.l++] = '!';
                }
                opts.ks_ins_seq.s [opts.ks_ins_seq.l] = 0;
                opts.ks_ins_qual.s[opts.ks_ins_qual.l] = 0;
            }
        }
        if (opts.last_tid >= 0)
            dump_fastq(&opts, sam_hdr_tid2name(opts.h, opts.last_tid),
                       opts.ks_ins_seq.s,  opts.ks_ins_seq.l,
                       opts.ks_ins_qual.s, opts.ks_ins_qual.l);
//        if (consensus_loop(&opts) < 0) {
//            print_error_errno("consensus", "Failed");
//            goto err;
//        }
    }

    ret = 0;

 err:
    if (opts.iter)
        hts_itr_destroy(opts.iter);
    if (opts.idx)
        hts_idx_destroy(opts.idx);

    if (opts.fp && sam_close(opts.fp) < 0) {
        print_error_errno("consensus", "Closing input file \"%s\"",
                          argv[optind]);
        ret = 1;
    }

    if (opts.h)
        sam_hdr_destroy(opts.h);
    sam_global_args_free(&ga);

    if (opts.fp_out && opts.fp_out != stdout)
        ret |= fclose(opts.fp_out) != 0;
    else
        ret |= fflush(stdout) != 0;

    ks_free(&opts.ks_line);
    ks_free(&opts.ks_ins_seq);
    ks_free(&opts.ks_ins_qual);

    if (ret)
        print_error("consensus", "failed");

    return ret;
}
