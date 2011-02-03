#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include "bam.h"

#define MAX_VARS 256

typedef struct {
	int8_t seq[MAX_VARS]; // TODO: change to dynamic memory allocation!
	int vpos, vlen;
} rseq_t, *rseq_p;

#define rseq_lt(a,b) ((a)->vpos < (b)->vpos)

#include "khash.h"
KHASH_MAP_INIT_INT64(64, rseq_t)

typedef khash_t(64) nseq_t;

#include "ksort.h"
KSORT_INIT(rseq, rseq_p, rseq_lt)

static int min_varQ = 40, min_mapQ = 10, var_len = 3;
static char nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

static inline uint64_t X31_hash_string(const char *s)
{
	uint64_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}

static void count_aux(int l, const uint8_t *seq, float *cnt)
{
	int i, j, n_ambi;
	uint32_t z, x;
	double y;
	for (i = n_ambi = 0; i < l; ++i)
		if (seq[i] == 0) ++n_ambi;
	if (n_ambi >= l - 1) return; // do nothing if too many ambiguous bases
	y = 1. / (1<<n_ambi);
	for (x = 0; x <= (1u<<n_ambi)-1; ++x) {
		for (i = j = 0, z = 0; i < l; ++i) {
			int c = seq[i]? seq[i] - 1 : (x>>j&1);
			z = z<<1 | c;
		}
		cnt[z] += y;
	}
}

static float **count_all(int l, int vpos, const nseq_t *hash)
{
	khint_t k;
	int i, j;
	uint8_t *seq;
	float **cnt;
	seq = calloc(l, 1);
	cnt = calloc(vpos, sizeof(void*));
	for (i = 0; i < vpos; ++i) cnt[i] = calloc(1<<l, sizeof(float));
	for (k = 0; k < kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			rseq_t *p = &kh_val(hash, k);
			if (p->vlen == 1) continue; // no phasing information
			for (j = 1; j < p->vlen; ++j) {
				for (i = 0; i < l; ++i)
					seq[i] = j < l - 1 - i? 0 : p->seq[j - (l - 1 - i)];
				count_aux(l, seq, cnt[p->vpos + j]);
			}
		}
	}
	free(seq);
	return cnt;
}

static void phase(int vpos, uint64_t *cns, nseq_t *hash)
{
	int i, j, n_seqs = kh_size(hash);
	khint_t k;
	rseq_t **seqs;
	float **cnt;
	cnt = count_all(var_len, vpos, hash);
	for (i = 0; i < vpos; ++i) {
		printf("%d", i);
		for (j = 0; j < 1<<var_len; ++j)
			printf("\t%.2f", cnt[i][j]);
		printf("\n");
	}
	for (i = 0; i < vpos; ++i) free(cnt[i]);
	free(cnt);
	seqs = calloc(n_seqs, sizeof(void*));
	for (k = 0, i = 0; k < kh_end(hash); ++k) 
		if (kh_exist(hash, k)) seqs[i++] = &kh_val(hash, k);
	ks_introsort_rseq(n_seqs, seqs);
	if (1) {
		for (i = 0; i < n_seqs; ++i) {
			printf("%d\t%d\t%d\t", (int)(cns[seqs[i]->vpos]>>32), seqs[i]->vpos, seqs[i]->vlen);
			for (j = 0; j < seqs[i]->vlen; ++j)
				putchar('0' + seqs[i]->seq[j]);
			putchar('\n');
		}
	}
//	for (i = 0; i < vpos; ++i) printf("%d\t%c\t%d\t%c\t%d\n", (int)(cns[i]>>32) + 1, "ACGT"[cns[i]&3], (int)(cns[i]&0xffff)>>2, "ACGT"[cns[i]>>16&3], (int)(cns[i]>>16&0xffff)>>2);
	free(seqs);
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
			phase(vpos, cns, seqs);
			lasttid = tid;
			vpos = 0;
		}
		cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
		for (i = 0; i < n; ++i) { // check if there are variants
			const bam_pileup1_t *p = plp + i;
			uint8_t *seq = bam1_seq(p->b);
			uint8_t *qual = bam1_qual(p->b);
			if (p->is_del || p->is_refskip) continue;
			if (p->b->core.qual < min_mapQ) continue;
			cnt[(int)nt16_nt4_table[(int)bam1_seqi(seq, p->qpos)]] += qual[p->qpos];
		}
		for (i = 0; i < 4; ++i) {
			if (cnt[i] >= 1<<14) cnt[i] = (1<<14) - 1;
			cnt[i] = cnt[i]<<2 | i; // cnt[i] is 16-bit at most
		}
		for (i = 1; i < 4; ++i) // insertion sort
			for (j = i; j > 0 && cnt[j] > cnt[j-1]; --j)
				tmp = cnt[j], cnt[j] = cnt[j-1], cnt[j-1] = tmp;
		if (cnt[1]>>2 <= min_varQ) continue; // not a variant
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
			rseq_t *r;
			if (p->is_del || p->is_refskip) continue;
			if (p->b->core.qual < min_mapQ) continue;
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
				if (vpos - r->vpos + 1 < MAX_VARS) {
					r->vlen = vpos - r->vpos + 1;
					r->seq[r->vlen-1] = c;
				}
			} else r->vpos = vpos, r->vlen = 1, r->seq[0] = c; // absent
		}
		++vpos;
	}
	phase(vpos, cns, seqs);
	bam_header_destroy(h);
	bam_plp_destroy(iter);
	bam_close(fp);
	kh_destroy(64, seqs);
	free(cns);
	return 0;
}
