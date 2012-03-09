/* To compile, copy this file to the samtools source code directory and compile with:
     gcc -g -O2 -Wall chk_indel_rg.c -o chk_indel_rg -Wall -I. -L. -lbam -lz
*/

#include <string.h>
#include "bam.h"

typedef struct {
	long cnt[4]; // short:ins, short:del, long:ins, long:del
} rgcnt_t;

#include "khash.h"
KHASH_MAP_INIT_STR(rgcnt, rgcnt_t)

#define MAX_LEN 127
#define Q_THRES 10
#define L_THRES 6 // short: <=L_THRES; otherwise long

int main(int argc, char *argv[])
{
	bamFile fp;
	bam1_t *b;
	int i, x;
	khash_t(rgcnt) *h;
	khint_t k;

	if (argc == 1) {
		fprintf(stderr, "Usage: chk_indel_rg <in.bam>\n\n");
		fprintf(stderr, "Output: filename, RG, #ins-in-short-homopolymer, #del-in-short, #ins-in-long, #del-in-long\n");
		return 1;
	}

	h = kh_init(rgcnt);
	fp = bam_open(argv[1], "r");
	bam_header_destroy(bam_header_read(fp)); // we do not need the header
	b = bam_init1();

	while (bam_read1(fp, b) >= 0) {
		if (b->core.n_cigar >= 3 && b->core.qual >= Q_THRES) {
			const uint8_t *seq;
			const uint32_t *cigar = bam1_cigar(b);
			char *rg;
			for (i = 0; i < b->core.n_cigar; ++i) // check if there are 1bp indels
				if (bam_cigar_oplen(cigar[i]) == 1 && (bam_cigar_op(cigar[i]) == BAM_CDEL || bam_cigar_op(cigar[i]) == BAM_CINS))
					break;
			if (i == b->core.n_cigar) continue; // no 1bp ins or del
			if ((rg = (char*)bam_aux_get(b, "RG")) == 0) continue; // no RG tag
			seq = bam1_seq(b);
			for (i = x = 0; i < b->core.n_cigar; ++i) {
				int op = bam_cigar_op(cigar[i]);
				if (bam_cigar_oplen(cigar[i]) == 1 && (op == BAM_CDEL || op == BAM_CINS)) {
					int c, j, hrun, which;
					c = bam1_seqi(seq, x);
					for (j = x + 1, hrun = 0; j < b->core.l_qseq; ++j, ++hrun) // calculate the hompolymer run length
						if (bam1_seqi(seq, j) != c) break;
					k = kh_get(rgcnt, h, rg + 1);
					if (k == kh_end(h)) { // absent
						char *key = strdup(rg + 1);
						k = kh_put(rgcnt, h, key, &c);
						memset(&kh_val(h, k), 0, sizeof(rgcnt_t));
					}
					which = (hrun <= L_THRES? 0 : 1)<<1 | (op == BAM_CINS? 0 : 1);
					++kh_val(h, k).cnt[which];
				}
				if (bam_cigar_type(op)&1) ++x;
			}
		}
	}

	for (k = 0; k != kh_end(h); ++k) {
		if (!kh_exist(h, k)) continue;
		printf("%s\t%s", argv[1], kh_key(h, k));
		for (i = 0; i < 4; ++i)
			printf("\t%ld", kh_val(h, k).cnt[i]);
		putchar('\n');
		free((char*)kh_key(h, k));
	}

	bam_destroy1(b);
	bam_close(fp);
	kh_destroy(rgcnt, h);
	return 0;
}
