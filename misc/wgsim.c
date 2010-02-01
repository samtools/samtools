/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

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

/* Contact: Heng Li <lh3@sanger.ac.uk> */

/* This program is separated from maq's read simulator with Colin
 * Hercus' modification to allow longer indels. Colin is the chief
 * developer of novoalign. */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>

#define PACKAGE_VERSION "0.2.3"

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

const int nst_color_space_table[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4, 4, 4};

/* Simple normal random number generator, copied from genran.c */

double ran_normal()
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

/* FASTA parser, copied from seq.c */

typedef struct {
	int l, m; /* length and maximum buffer size */
	unsigned char *s; /* sequence */
} seq_t;

#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).m = 0

static int SEQ_BLOCK_SIZE = 512;

void seq_set_block_size(int size)
{
	SEQ_BLOCK_SIZE = size;
}

int seq_read_fasta(FILE *fp, seq_t *seq, char *locus, char *comment)
{
	int c, l, max;
	char *p;
	
	c = 0;
	while (!feof(fp) && fgetc(fp) != '>');
	if (feof(fp)) return -1;
	p = locus;
	while (!feof(fp) && (c = fgetc(fp)) != ' ' && c != '\t' && c != '\n')
		if (c != '\r') *p++ = c;
	*p = '\0';
	if (comment) {
		p = comment;
		if (c != '\n') {
			while (!feof(fp) && ((c = fgetc(fp)) == ' ' || c == '\t'));
			if (c != '\n') {
				*p++ = c;
				while (!feof(fp) && (c = fgetc(fp)) != '\n')
					if (c != '\r') *p++ = c;
			}
		}
		*p = '\0';
	} else if (c != '\n') while (!feof(fp) && fgetc(fp) != '\n');
	l = 0; max = seq->m;
	while (!feof(fp) && (c = fgetc(fp)) != '>') {
		if (isalpha(c) || c == '-' || c == '.') {
			if (l + 1 >= max) {
				max += SEQ_BLOCK_SIZE;
				seq->s = (unsigned char*)realloc(seq->s, sizeof(char) * max);
			}
			seq->s[l++] = (unsigned char)c;
		}
	}
	if (c == '>') ungetc(c,fp);
	seq->s[l] = 0;
	seq->m = max; seq->l = l;
	return l;
}

/* Error-checking open, copied from utils.c */

#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r"))? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	return fp;
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
static double INDEL_FRAC = 0.1;
static double INDEL_EXTEND = 0.3;
static int IS_SOLID = 0;
static int SHOW_MM_INFO = 1;

void maq_mut_diref(const seq_t *seq, int is_hap, mutseq_t *hap1, mutseq_t *hap2)
{
	int i, deleting = 0;
	mutseq_t *ret[2];

	ret[0] = hap1; ret[1] = hap2;
	ret[0]->l = seq->l; ret[1]->l = seq->l;
	ret[0]->m = seq->m; ret[1]->m = seq->m;
	ret[0]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
	ret[1]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
	for (i = 0; i != seq->l; ++i) {
		int c;
		c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
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
void maq_print_mutref(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2)
{
	int i;
	for (i = 0; i != seq->l; ++i) {
		int c[3];
		c[0] = nst_nt4_table[(int)seq->s[i]];
		c[1] = hap1->s[i]; c[2] = hap2->s[i];
		if (c[0] >= 4) continue;
		if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
			printf("%s\t%d\t", name, i+1);
			if (c[1] == c[2]) { // hom
				if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
					printf("%c\t%c\t-\n", "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
				} else if ((c[1]&mutmsk) == DELETE) { // del
					printf("%c\t-\t-\n", "ACGTN"[c[0]]);
				} else if (((c[1] & mutmsk) >> 12) <= 5) { // ins
					printf("-\t");
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while(n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        n--;
                    }
                    printf("\t-\n");
				}  else assert(0);
			} else { // het
				if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
					printf("%c\t%c\t+\n", "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
				} else if ((c[1]&mutmsk) == DELETE) {
					printf("%c\t-\t+\n", "ACGTN"[c[0]]);
				} else if ((c[2]&mutmsk) == DELETE) {
					printf("%c\t-\t+\n", "ACGTN"[c[0]]);
				} else if (((c[1] & mutmsk) >> 12) <= 4) { // ins1
					printf("-\t");
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        n--;
                    }
                    printf("\t+\n");
				} else if (((c[2] & mutmsk) >> 12) <= 5) { // ins2
					printf("-\t");
                    int n = (c[2]&mutmsk) >> 12, ins = c[2] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
				} else assert(0);
			}
		}
	}
}

void wgsim_core(FILE *fpout1, FILE *fpout2, FILE *fp_fa, int is_hap, uint64_t N, int dist, int std_dev, int size_l, int size_r)
{
	seq_t seq;
    mutseq_t rseq[2];
	uint64_t tot_len, ii;
	int i, l, n_ref;
	char name[256], *qstr;
	int size[2], Q;
	uint8_t *tmp_seq[2];
    mut_t *target;

	INIT_SEQ(seq);
	srand48(time(0));
	seq_set_block_size(0x1000000);
	l = size_l > size_r? size_l : size_r;
	qstr = (char*)calloc(l+1, 1);
	tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
	tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
	size[0] = size_l; size[1] = size_r;

	Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

	tot_len = n_ref = 0;
	while ((l = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
		tot_len += l;
		++n_ref;
	}
	fprintf(stderr, "[wgsim_core] %d sequences, total length: %llu\n", n_ref, (long long)tot_len);
	rewind(fp_fa);

	while ((l = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
		uint64_t n_pairs = (uint64_t)((long double)l / tot_len * N + 0.5);
		if (l < dist + 3 * std_dev) {
			fprintf(stderr, "[wgsim_core] kkip sequence '%s' as it is shorter than %d!\n", name, dist + 3 * std_dev);
			continue;
		}

		// generate mutations and print them out
		maq_mut_diref(&seq, is_hap, rseq, rseq+1);
		maq_print_mutref(name, &seq, rseq, rseq+1);

		for (ii = 0; ii != n_pairs; ++ii) { // the core loop
			double ran;
			int d, pos, s[2], is_flip = 0;
			int n_sub[2], n_indel[2], n_err[2], ext_coor[2], j, k;
			FILE *fpo[2];

			do { // avoid boundary failure
				ran = ran_normal();
				ran = ran * std_dev + dist;
				d = (int)(ran + 0.5);
				pos = (int)((l - d + 1) * drand48());
			} while (pos < 0 || pos >= seq.l || pos + d - 1 >= seq.l);

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

#define __gen_read(x, start, iter) do {									\
				for (i = (start), k = 0, ext_coor[x] = -10; i >= 0 && i < seq.l && k < s[x]; iter) {	\
					int c = target[i], mut_type = c & mutmsk;			\
					if (ext_coor[x] < 0) {								\
						if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
						ext_coor[x] = i;								\
					}													\
					if (mut_type == DELETE) ++n_indel[x];				\
					else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
						tmp_seq[x][k++] = c & 0xf;						\
						if (mut_type == SUBSTITUTE) ++n_sub[x];			\
					} else {											\
						int n, ins;										\
						++n_indel[x];									\
						tmp_seq[x][k++] = c & 0xf;						\
						for (n = mut_type>>12, ins = c>>4; n > 0 && k < s[x]; --n, ins >>= 2) \
							tmp_seq[x][k++] = ins & 0x3;				\
					}													\
				}														\
				if (k != s[x]) ext_coor[x] = -10;						\
			} while (0)

			if (!IS_SOLID) {
				__gen_read(0, pos, ++i);
				__gen_read(1, pos + d - 1, --i);
				for (k = 0; k < s[1]; ++k) tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4; // complement
			} else {
				int c1, c2, c;
				++s[0]; ++s[1]; // temporarily increase read length by 1
				if (is_flip) { // RR pair
					__gen_read(0, pos + s[0], --i);
					__gen_read(1, pos + d - 1, --i);
				} else { // FF pair
					__gen_read(0, pos, ++i);
					__gen_read(1, pos + d - 1 - s[1], ++i);
					++ext_coor[0]; ++ext_coor[1];
				}
				// change to color sequence: (0,1,2,3) -> (A,C,G,T)
				for (j = 0; j < 2; ++j) {
					c1 = tmp_seq[j][0];
					for (i = 1; i < s[j]; ++i) {
						c2 = tmp_seq[j][i];
						c = (c1 >= 4 || c2 >= 4)? 4 : nst_color_space_table[(1<<c1)|(1<<c2)];
						tmp_seq[j][i-1] = c;
						c1 = c2;
					}
				}
				--s[0]; --s[1]; // change back
			}
			if (ext_coor[0] < 0 || ext_coor[1] < 0) { // fail to generate the read(s)
				--ii;
				continue;
			}

			// generate sequencing errors
			for (j = 0; j < 2; ++j) {
				for (i = 0; i < s[j]; ++i) {
					int c = tmp_seq[j][i];
					if (c >= 4) c = 4; // actually c should be never larger than 4 if everything is correct
					else if (drand48() < ERR_RATE) {
						c = (c + (int)(drand48() * 3.0 + 1)) & 3;
						++n_err[j];
					}
					tmp_seq[j][i] = c;
				}
			}

			// print
			for (j = 0; j < 2; ++j) {
				for (i = 0; i < s[j]; ++i) qstr[i] = Q;
				qstr[i] = 0;
				if (SHOW_MM_INFO) {
					fprintf(fpo[j], "@%s_%u_%u_%d:%d:%d_%d:%d:%d_%llx/%d\n", name, ext_coor[0]+1, ext_coor[1]+1,
							n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
							(long long)ii, j==0? is_flip+1 : 2-is_flip);
				} else {
					fprintf(fpo[j], "@%s_%u_%u_%llx/%d %d:%d:%d_%d:%d:%d\n", name, ext_coor[0]+1, ext_coor[1]+1,
							(long long)ii, j==0? is_flip+1 : 2-is_flip,
							n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1]);
				}
				for (i = 0; i < s[j]; ++i)
					fputc("ACGTN"[(int)tmp_seq[j][i]], fpo[j]);
				fprintf(fpo[j], "\n+\n%s\n", qstr);
			}
		}
		free(rseq[0].s); free(rseq[1].s);
	}
	free(seq.s); free(qstr);
	free(tmp_seq[0]); free(tmp_seq[1]);
}

static int simu_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: wgsim (short read simulator)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
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
	fprintf(stderr, "         -c            generate reads in color space (SOLiD reads)\n");
	fprintf(stderr, "         -C            show mismatch info in comment rather than read name\n");
	fprintf(stderr, "         -h            haplotype mode\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Note: For SOLiD reads, the first read is F3 and the second is R3.\n\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int64_t N;
	int dist, std_dev, c, size_l, size_r, is_hap = 0;
	FILE *fpout1, *fpout2, *fp_fa;

	N = 1000000; dist = 500; std_dev = 50;
	size_l = size_r = 70;
	while ((c = getopt(argc, argv, "e:d:s:N:1:2:r:R:hX:cC")) >= 0) {
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
		case 'c': IS_SOLID = 1; break;
		case 'C': SHOW_MM_INFO = 0; break;
		case 'h': is_hap = 1; break;
		}
	}
	if (argc - optind < 3) return simu_usage();
	fp_fa = (strcmp(argv[optind+0], "-") == 0)? stdin : xopen(argv[optind+0], "r");
	fpout1 = xopen(argv[optind+1], "w");
	fpout2 = xopen(argv[optind+2], "w");
	wgsim_core(fpout1, fpout2, fp_fa, is_hap, N, dist, std_dev, size_l, size_r);

	fclose(fpout1); fclose(fpout2); fclose(fp_fa);
	return 0;
}
