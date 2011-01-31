#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include "bam.h"

typedef struct {
	uint64_t seq;
	int vpos; // 4 bytes wasted on 64-bit machines
} rseq_t;

#include "khash.h"
KHASH_MAP_INIT_INT64(64, rseq_t)

typedef khash_t(64) nseq_t;

static int min_Q = 20;
static char nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

static inline uint64_t X31_hash_string(const char *s)
{
	uint64_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}

static void phase(int vpos, uint64_t *cns)
{
	int i;
	for (i = 0; i < vpos; ++i) {
		printf("%d\t%c\t%d\t%c\t%d\n", (int)(cns[i]>>32) + 1, "ACGT"[cns[i]&3], (int)(cns[i]&0xffff)>>2, "ACGT"[cns[i]>>16&3], (int)(cns[i]>>16&0xffff)>>2);
	}
}

int main_phase(int argc, char *argv[])
{
	bamFile fp;
	int c, tid, pos, vpos = 0, n, lasttid = -1, max_vpos = 0;
	const bam_pileup1_t *plp;
	bam_plp_t iter;
	bam_header_t *h;
	nseq_t *seqs;
	uint64_t *cns = 0;

	while ((c = getopt(argc, argv, "")) >= 0) {
		switch (c) {
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samtools phase <in.bam>\n");
		return 1;
	}
	fp = bam_open(argv[optind], "r");
	h = bam_header_read(fp);
	iter = bam_plp_init((bam_plp_auto_f)bam_read1, fp);

	seqs = kh_init(64);
	while ((plp = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {
		int i, j, c, cnt[4], tmp;
		if (tid < 0) break;
		if (tid != lasttid) { // change of chromosome
			phase(vpos, cns);
			lasttid = tid;
			vpos = 0;
		}
		cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
		for (i = 0; i < n; ++i) { // check if there are variants
			const bam_pileup1_t *p = plp + i;
			uint8_t *seq = bam1_seq(p->b);
			uint8_t *qual = bam1_qual(p->b);
			if (p->is_del || p->is_refskip) continue;
			if (p->b->core.qual < min_Q || qual[p->qpos] < min_Q) continue;
			++cnt[(int)nt16_nt4_table[(int)bam1_seqi(seq, p->qpos)]];
		}
		for (i = 0; i < 4; ++i) {
			if (cnt[i] >= 1<<14) cnt[i] = (1<<14) - 1;
			cnt[i] = cnt[i]<<2 | i; // cnt[i] is 16-bit at most
		}
		for (i = 1; i < 4; ++i) // insertion sort
			for (j = i; j > 0 && cnt[j] > cnt[j-1]; --j)
				tmp = cnt[j], cnt[j] = cnt[j-1], cnt[j-1] = tmp;
		if (cnt[1]>>2 == 0) continue; // not a variant
		if (vpos == max_vpos) {
			max_vpos = max_vpos? max_vpos<<1 : 128;
			cns = realloc(cns, max_vpos * 8);
		}
		cns[vpos] = (uint64_t)pos<<32 | cnt[1] << 16 | cnt[0];
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = plp + i;
			uint64_t key;
			khint_t k;
			uint8_t *seq = bam1_seq(p->b);
			uint8_t *qual = bam1_qual(p->b);
			rseq_t *r;
			if (p->is_del || p->is_refskip) continue;
			if (p->b->core.qual < min_Q || qual[p->qpos] < min_Q) continue;
			// get the base code
			c = nt16_nt4_table[(int)bam1_seqi(seq, p->qpos)];
			if (c > 3) c = 0;
			else if (c == (cnt[0]&3)) c = 1;
			else if (c == (cnt[1]&3)) c = 2;
			else c = 0; // TODO: perhaps c=3 is better? Watch out other dependencies!
			// write to seqs
			key = X31_hash_string(bam1_qname(p->b));
			k = kh_put(64, seqs, key, &tmp);
			r = &kh_val(seqs, k);
			if (tmp == 0) { // present in the hash table
				if (vpos - r->vpos < 32) // 32 bases at the maximum
					r->seq |= c << (vpos - r->vpos)*2;
			} else r->vpos = vpos, r->seq = c; // absent
		}
		++vpos;
	}
	phase(vpos, cns);
	bam_header_destroy(h);
	bam_plp_destroy(iter);
	bam_close(fp);
	kh_destroy(64, seqs);
	free(cns);
	return 0;
}
