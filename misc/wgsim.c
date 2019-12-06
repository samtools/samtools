/* The MIT License

   Copyright (c) 2008, 2018 Genome Research Ltd (GRL).
                 2011 Heng Li <lh3@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* This program is separated from maq's read simulator with Colin
 * Hercus' modification to allow longer indels. */

#include <config.h>

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include "../version.h"
#include "htslib/kseq.h"
#include "htslib/hts_os.h"
KSEQ_INIT(gzFile, gzread)

const uint8_t nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/* Simple normal random number generator, copied from genran.c */

double ran_normal(void)
{
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;
    if (iset == 0) {
        do {
            v1 = 2.0 * drand48() - 1.0;
            v2 = 2.0 * drand48() - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }
}

/* wgsim */

enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;

static double ERR_RATE = 0.02;
static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.15;
static double INDEL_EXTEND = 0.3;
static double MAX_N_RATIO = 0.05;

void wgsim_mut_diref(const kseq_t *ks, int is_hap, mutseq_t *hap1, mutseq_t *hap2)
{
    int i, deleting = 0;
    mutseq_t *ret[2];

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    for (i = 0; i != ks->seq.l; ++i) {
        int c;
        c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (deleting) {
            if (drand48() < INDEL_EXTEND) {
                if (deleting & 1) ret[0]->s[i] |= DELETE;
                if (deleting & 2) ret[1]->s[i] |= DELETE;
                continue;
            } else deleting = 0;
        }
        if (c < 4 && drand48() < MUT_RATE) { // mutation
            if (drand48() >= INDEL_FRAC) { // substitution
                double r = drand48();
                c = (c + (int)(r * 3.0 + 1)) & 3;
                if (is_hap || drand48() < 0.333333) { // hom
                    ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                } else { // het
                    ret[drand48()<0.5?0:1]->s[i] = SUBSTITUTE|c;
                }
            } else { // indel
                if (drand48() < 0.5) { // deletion
                    if (is_hap || drand48() < 0.333333) { // hom-del
                        ret[0]->s[i] = ret[1]->s[i] = DELETE;
                        deleting = 3;
                    } else { // het-del
                        deleting = drand48()<0.5?1:2;
                        ret[deleting-1]->s[i] = DELETE;
                    }
                } else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        num_ins++;
                        ins = (ins << 2) | (int)(drand48() * 4.0);
                    } while (num_ins < 4 && drand48() < INDEL_EXTEND);

                    if (is_hap || drand48() < 0.333333) { // hom-ins
                        ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    } else { // het-ins
                        ret[drand48()<0.5?0:1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    }
                }
            }
        }
    }
}
void wgsim_print_mutref(const char *name, const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2)
{
    int i, j = 0; // j keeps the end of the last deletion
    for (i = 0; i != ks->seq.l; ++i) {
        int c[3];
        c[0] = nst_nt4_table[(int)ks->seq.s[i]];
        c[1] = hap1->s[i]; c[2] = hap2->s[i];
        if (c[0] >= 4) continue;
        if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
            if (c[1] == c[2]) { // hom
                if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t-\n", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
                } else if ((c[1]&mutmsk) == DELETE) { // del
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+1);
                        for (j = i; j < ks->seq.l && hap1->s[j] == hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4) { // ins
                    printf("%s\t%d\t-\t", name, i+1);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t-\n");
                } // else: deleted base in a long deletion
            } else { // het
                if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t+\n", name, i+1, "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
                } else if ((c[1]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+1);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if ((c[2]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+1);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap2->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4 && ((c[1] & mutmsk) >> 12) > 0) { // ins1
                    printf("%s\t%d\t-\t", name, i+1);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
                } else if (((c[2] & mutmsk) >> 12) <= 4 || ((c[2] & mutmsk) >> 12) > 0) { // ins2
                    printf("%s\t%d\t-\t", name, i+1);
                    int n = (c[2]&mutmsk) >> 12, ins = c[2] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
                } // else: deleted base in a long deletion
            }
        }
    }
}

void wgsim_core(FILE *fpout1, FILE *fpout2, const char *fn, int is_hap, uint64_t N, int dist, int std_dev, int size_l, int size_r)
{
    kseq_t *ks;
    mutseq_t rseq[2];
    gzFile fp_fa;
    uint64_t tot_len, ii;
    int i, l, n_ref;
    char *qstr;
    int size[2], Q, max_size;
    uint8_t *tmp_seq[2];
    mut_t *target;
    int max_loop, max_loop_err = 0;

    fp_fa = gzopen(fn, "r");
    if (fp_fa == NULL) {
        fprintf(stderr, "[wgsim] can't open '%s': %s\n", fn, strerror(errno));
        return;
    }

    l = size_l > size_r? size_l : size_r;
    qstr = (char*)calloc(l+1, 1);
    tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
    tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
    size[0] = size_l; size[1] = size_r;
    max_size = size_l > size_r? size_l : size_r;

    Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

    ks = kseq_init(fp_fa);
    tot_len = n_ref = 0;
    fprintf(stderr, "[%s] calculating the total length of the reference sequence...\n", __func__);
    while ((l = kseq_read(ks)) >= 0) {
        tot_len += l;
        ++n_ref;
    }
    fprintf(stderr, "[%s] %d sequences, total length: %llu\n", __func__, n_ref, (long long)tot_len);
    kseq_destroy(ks);
    gzclose(fp_fa);

    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);
    while ((l = kseq_read(ks)) >= 0) {
        uint64_t n_pairs = (uint64_t)((long double)l / tot_len * N + 0.5);
        if (l < dist + 3 * std_dev) {
            fprintf(stderr, "[%s] skip sequence '%s' as it is shorter than %d!\n", __func__, ks->name.s, dist + 3 * std_dev);
            continue;
        }

        // generate mutations and print them out
        wgsim_mut_diref(ks, is_hap, rseq, rseq+1);
        wgsim_print_mutref(ks->name.s, ks, rseq, rseq+1);

        for (ii = 0; max_loop = 1000, ii != n_pairs; ++ii) { // the core loop
            double ran;
            int d, pos, s[2], is_flip;
            int n_sub[2], n_indel[2], n_err[2], ext_coor[2], j, k;
            FILE *fpo[2];

        try_again:
            is_flip = 0;

            do { // avoid boundary failure
                ran = ran_normal();
                ran = ran * std_dev + dist;
                d = (int)(ran + 0.5);
                d = d > max_size? d : max_size;
                pos = (int)((l - d + 1) * drand48());
            } while (pos < 0 || pos >= ks->seq.l || pos + d - 1 >= ks->seq.l);

            // flip or not
            if (drand48() < 0.5) {
                fpo[0] = fpout1; fpo[1] = fpout2;
                s[0] = size[0]; s[1] = size[1];
            } else {
                fpo[1] = fpout1; fpo[0] = fpout2;
                s[1] = size[0]; s[0] = size[1];
                is_flip = 1;
            }

            // generate the read sequences
            target = rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
            n_sub[0] = n_sub[1] = n_indel[0] = n_indel[1] = n_err[0] = n_err[1] = 0;

            // forward read
            for (i = pos, k = 0, ext_coor[0] = -10; i >= 0 && i < ks->seq.l && k < s[0]; ++i) {
                int c = target[i], mut_type = c & mutmsk;
                if (ext_coor[0] < 0) {
                    if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue;
                    ext_coor[0] = i;
                }
                if (mut_type == DELETE) ++n_indel[0];
                else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) {
                    tmp_seq[0][k++] = c & 0xf;
                    if (mut_type == SUBSTITUTE) ++n_sub[0];
                } else {
                    int n, ins;
                    ++n_indel[0];
                    tmp_seq[0][k++] = c & 0xf;
                    for (n = mut_type>>12, ins = c>>4; n > 0 && k < s[0]; --n, ins >>= 2)
                        tmp_seq[0][k++] = ins & 0x3;
                }
            }
            if (k != s[0]) ext_coor[0] = -10;

            // reverse read
            for (i = pos + d - 1, k = 0, ext_coor[1] = -10; i >= 0 && i < ks->seq.l && k < s[1]; --i) {
                int c = target[i], mut_type = c & mutmsk;
                if (ext_coor[1] < 0) {
                    if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue;
                    ext_coor[1] = i;
                }
                if (mut_type == DELETE) ++n_indel[1];
                else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) {
                    tmp_seq[1][k++] = c & 0xf;
                    if (mut_type == SUBSTITUTE) ++n_sub[1];
                } else {
                    int n, ins;
                    ++n_indel[1];
                    for (n = mut_type>>12, ins = c>>4; n > 0 && k < s[1];)
                        tmp_seq[1][k++] = (ins>>(2*--n)) & 0x3;
                    tmp_seq[1][k++] = c & 0xf;
                }
            }
            if (k != s[1]) ext_coor[1] = -10;



            for (k = 0; k < s[1]; ++k) tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4; // complement
            if (ext_coor[0] < 0 || ext_coor[1] < 0) { // fail to generate the read(s)
                --ii;
                continue;
            }

            // generate sequencing errors
            for (j = 0; j < 2; ++j) {
                int n_n = 0;
                for (i = 0; i < s[j]; ++i) {
                    int c = tmp_seq[j][i];
                    if (c >= 4) { // actually c should be never larger than 4 if everything is correct
                        c = 4;
                        ++n_n;
                    } else if (drand48() < ERR_RATE) {
                        // c = (c + (int)(drand48() * 3.0 + 1)) & 3; // random sequencing errors
                        c = (c + 1) & 3; // recurrent sequencing errors
                        ++n_err[j];
                    }
                    tmp_seq[j][i] = c;
                }
                if ((double)n_n / s[j] > MAX_N_RATIO) break;
            }
            if (j < 2) { // too many ambiguous bases on one of the reads
                if (max_loop--) goto try_again;
                if (!max_loop_err) {
                    fprintf(stderr, "Failed to produce a sequence with insufficient Ns. "
                            "Omitting some sequence-pairs\n");
                    max_loop_err = 1;
                }
                continue;
            }

            // print
            for (j = 0; j < 2; ++j) {
                for (i = 0; i < s[j]; ++i) qstr[i] = Q;
                qstr[i] = 0;
                fprintf(fpo[j], "@%s_%u_%u_%d:%d:%d_%d:%d:%d_%llx/%d\n", ks->name.s, ext_coor[0]+1, ext_coor[1]+1,
                        n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
                        (long long)ii, j==0? is_flip+1 : 2-is_flip);
                for (i = 0; i < s[j]; ++i)
                    fputc("ACGTN"[(int)tmp_seq[j][i]], fpo[j]);
                fprintf(fpo[j], "\n+\n%s\n", qstr);
            }
        }
        free(rseq[0].s); free(rseq[1].s);
    }
    kseq_destroy(ks);
    gzclose(fp_fa);
    free(qstr);
    free(tmp_seq[0]); free(tmp_seq[1]);
}

static int simu_usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: wgsim (short read simulator)\n");
    fprintf(stderr, "Version: %s\n", SAMTOOLS_VERSION);
    fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
    fprintf(stderr, "Usage:   wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>\n\n");
    fprintf(stderr, "Options: -e FLOAT      base error rate [%.3f]\n", ERR_RATE);
    fprintf(stderr, "         -d INT        outer distance between the two ends [500]\n");
    fprintf(stderr, "         -s INT        standard deviation [50]\n");
    fprintf(stderr, "         -N INT        number of read pairs [1000000]\n");
    fprintf(stderr, "         -1 INT        length of the first read [70]\n");
    fprintf(stderr, "         -2 INT        length of the second read [70]\n");
    fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", MUT_RATE);
    fprintf(stderr, "         -R FLOAT      fraction of indels [%.2f]\n", INDEL_FRAC);
    fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", INDEL_EXTEND);
    fprintf(stderr, "         -S INT        seed for random generator [0, use the current time]\n");
    fprintf(stderr, "         -A FLOAT      discard if the fraction of ambiguous bases higher than FLOAT [%.2f]\n", MAX_N_RATIO);
    fprintf(stderr, "         -h            haplotype mode\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{
    int64_t N;
    int dist, std_dev, c, size_l, size_r, is_hap = 0;
    FILE *fpout1, *fpout2;
    int seed = 0;

    N = 1000000; dist = 500; std_dev = 50;
    size_l = size_r = 70;
    while ((c = getopt(argc, argv, "e:d:s:N:1:2:r:R:hX:S:A:")) >= 0) {
        switch (c) {
        case 'd': dist = atoi(optarg); break;
        case 's': std_dev = atoi(optarg); break;
        case 'N': N = atoi(optarg); break;
        case '1': size_l = atoi(optarg); break;
        case '2': size_r = atoi(optarg); break;
        case 'e': ERR_RATE = atof(optarg); break;
        case 'r': MUT_RATE = atof(optarg); break;
        case 'R': INDEL_FRAC = atof(optarg); break;
        case 'X': INDEL_EXTEND = atof(optarg); break;
        case 'A': MAX_N_RATIO = atof(optarg); break;
        case 'S': seed = atoi(optarg); break;
        case 'h': is_hap = 1; break;
        }
    }
    if (argc - optind < 3) return simu_usage();
    fpout1 = fopen(argv[optind+1], "w");
    fpout2 = fopen(argv[optind+2], "w");
    if (!fpout1 || !fpout2) {
        fprintf(stderr, "[wgsim] file open error\n");
        return 1;
    }
    if (seed <= 0) seed = time(0)&0x7fffffff;
    fprintf(stderr, "[wgsim] seed = %d\n", seed);
    srand48(seed);
    wgsim_core(fpout1, fpout2, argv[optind], is_hap, N, dist, std_dev, size_l, size_r);

    fclose(fpout1); fclose(fpout2);
    return 0;
}
