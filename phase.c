#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include "bam.h"

#define MAX_VARS 256

typedef struct {
	int8_t seq[MAX_VARS]; // TODO: change to dynamic memory allocation!
	int vpos, beg, end;
	uint32_t vlen:30, phase:2;
} rseq_t, *rseq_p;

#define rseq_lt(a,b) ((a)->beg < (b)->beg)

#include "khash.h"
KHASH_MAP_INIT_INT64(64, rseq_t)

typedef khash_t(64) nseq_t;

#include "ksort.h"
KSORT_INIT(rseq, rseq_p, rseq_lt)

static int min_varQ = 40, min_mapQ = 10, var_len = 5;
static int g_vpos_shift = 0;
static char nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

static inline uint64_t X31_hash_string(const char *s)
{
	uint64_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}

static void count1(int l, const uint8_t *seq, int *cnt)
{
	int i, j, n_ambi;
	uint32_t z, x;
	if (seq[l-1] == 0) return; // do nothing is the last base is ambiguous
	for (i = n_ambi = 0; i < l; ++i) // collect ambiguous bases
		if (seq[i] == 0) ++n_ambi;
	if (l - n_ambi <= 1) return; // only one SNP
	for (x = 0; x < 1u<<n_ambi; ++x) { // count
		for (i = j = 0, z = 0; i < l; ++i) {
			int c;
			if (seq[i]) c = seq[i] - 1;
			else {
				c = x>>j&1;
				++j;
			}
			z = z<<1 | c;
		}
		++cnt[z];
	}
}

static int **count_all(int l, int vpos, const nseq_t *hash)
{
	khint_t k;
	int i, j, **cnt;
	uint8_t *seq;
	seq = calloc(l, 1);
	cnt = calloc(vpos, sizeof(void*));
	for (i = 0; i < vpos; ++i) cnt[i] = calloc(1<<l, sizeof(int));
	for (k = 0; k < kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			rseq_t *p = &kh_val(hash, k);
			if (p->vlen == 1 || p->vpos >= vpos) continue; // no phasing information or out of region
			for (j = 1; j < p->vlen; ++j) {
				for (i = 0; i < l; ++i)
					seq[i] = j < l - 1 - i? 0 : p->seq[j - (l - 1 - i)];
				count1(l, seq, cnt[p->vpos + j]);
			}
		}
	}
	free(seq);
	return cnt;
}

static int8_t *dynaprog(int l, int vpos, int **w)
{
	int *f[2], *curr, *prev, max, i;
	int8_t **b, *h = 0;
	uint32_t x, z = 1u<<(l-1), mask = (1u<<l) - 1;
	f[0] = calloc(z, sizeof(int));
	f[1] = calloc(z, sizeof(int));
	b = calloc(vpos, sizeof(void*));
	prev = f[0]; curr = f[1];
	// fill the backtrack matrix
	for (i = 0; i < vpos; ++i) {
		int *wi = w[i], *tmp;
		int8_t *bi;
		bi = b[i] = calloc(z, 1);
		/* In the following, x is the current state, which is the
		 * lexicographically smaller local haplotype. xc is the complement of
		 * x, or the larger local haplotype; y0 and y1 are the two predecessors
		 * of x. */
		for (x = 0; x < z; ++x) { // x0 is the smaller 
			uint32_t y0, y1, xc;
			int c0, c1;
			xc = ~x&mask; y0 = x>>1; y1 = xc>>1;
			c0 = prev[y0] + wi[x] + wi[xc];
			c1 = prev[y1] + wi[x] + wi[xc];
			if (c0 > c1) bi[x] = 0, curr[x] = c0;
			else bi[x] = 1, curr[x] = c1;
		}
		tmp = prev; prev = curr; curr = tmp; // swap
	}
	{ // backtrack
		uint32_t max_x = 0;
		int which = 0;
		h = calloc(vpos, 1);
		for (x = 0, max = 0, max_x = 0; x < z; ++x)
			if (prev[x] > max) max = prev[x], max_x = x;
		for (i = vpos - 1, x = max_x; i >= 0; --i) {
			h[i] = which? (~x&1) : (x&1);
			which = b[i][x]? !which : which;
			x = b[i][x]? (~x&mask)>>1 : x>>1;
		}
	}
	// free
	for (i = 0; i < vpos; ++i) free(b[i]);
	free(f[0]); free(f[1]); free(b);
	return h;
}

static int filter(int vpos, int *const* cnt, const int8_t *path, uint64_t *cns, nseq_t *hash)
{
	int8_t *flt;
	int i, k, *map;
	khint_t j;
	uint32_t x0, x1, mask = (1<<var_len) - 1;
	flt = calloc(vpos, 1);
	map = calloc(vpos, sizeof(int));
	// get the list of sites to be filtered
	for (i = 1, x0 = x1 = 0; i < vpos; ++i) {
		int *ci = cnt[i];
		x0 = (x0<<1 | path[i]) & mask; x1 = ~x0 & mask;
		if (ci[x0] == 0 || ci[x1] == 0) flt[i] = 1; // no supporting fragment for either haplotype
	}
	// generate map[]
	for (i = k = 0; i < vpos; ++i) {
		if (flt[i]) map[i] = -1;
		else map[i] = k++;
	}
	// filter hash
	for (j = 0; j < kh_end(hash); ++j) {
		if (kh_exist(hash, j)) {
			rseq_t *s = &kh_val(hash, j);
			if (s->vpos >= vpos) continue;
			for (i = k = 0; i < s->vlen; ++i)
				if (flt[s->vpos + i] == 0 && i != k)
					s->seq[k++] = s->seq[i];
			if (k < 2) kh_del(64, hash, j);
			s->vlen = k; s->vpos = map[s->vpos];
		}
	}
	// filter cns
	for (i = k = 0; i < vpos; ++i)
		if (flt[i] == 0) cns[k++] = cns[i];
	free(flt); free(map);
	return k;
}

static void phase(const char *chr, int vpos, uint64_t *cns, nseq_t *hash)
{
	int **cnt, i, j, n_seqs = kh_size(hash), ori_vpos = vpos;
	khint_t k;
	rseq_t **seqs;
	int8_t *path = 0;
	if (vpos == 0) return;
	printf("BL\t%s\t%d\t%d\n", chr, (int)(cns[0]>>32), (int)(cns[vpos-1]>>32));
	//filter(vpos, cnt, cns, hash);
	cnt = count_all(var_len, vpos, hash);
	path = dynaprog(var_len, vpos, cnt);
	{
		uint32_t x0, x1, mask = (1<<var_len) - 1;
		for (i = 0, x0 = x1 = 0; i < vpos; ++i) {
			int8_t c[2];
			c[0] = cns[i]&3; c[1] = cns[i]>>16&3;
			x0 = (x0<<1 | path[i]) & mask ; x1 = ~x0 & mask;
			printf("VL\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\n", (int)(cns[i]>>32) + 1, i + g_vpos_shift + 1, "ACGT"[c[path[i]]], "ACGT"[c[1-path[i]]], path[i],
				cnt[i][x0], cnt[i][x0^1], cnt[i][x1], cnt[i][x1^1]);
			//for (j = 0; j < 1<<var_len; ++j) printf("%c%d", (j&1)? ' ' : '\t', cnt[i][j]);
		}
	}
	for (i = 0; i < vpos; ++i) free(cnt[i]);
	free(cnt); free(path);
	seqs = calloc(n_seqs, sizeof(void*));
	for (k = 0, i = 0; k < kh_end(hash); ++k) 
		if (kh_exist(hash, k) && kh_val(hash, k).vpos < vpos) seqs[i++] = &kh_val(hash, k);
	n_seqs = i;
	ks_introsort_rseq(n_seqs, seqs);
	for (i = 0; i < n_seqs; ++i) {
		rseq_t *s = seqs[i];
		printf("EV\t0\t%s\t%d\t40\t%dM\t*\t0\t0\t", chr, s->vpos + 1 + g_vpos_shift, s->vlen);
		for (j = 0; j < s->vlen; ++j) {
			uint32_t c = cns[s->vpos + j];
			if (s->seq[j] == 0) putchar('N');
			else putchar("ACGT"[s->seq[j] == 1? (c&3) : (c>>16&3)]);
		}
		printf("\t*\n");
	}
//	for (i = 0; i < vpos; ++i) printf("%d\t%c\t%d\t%c\t%d\n", (int)(cns[i]>>32) + 1, "ACGT"[cns[i]&3], (int)(cns[i]&0xffff)>>2, "ACGT"[cns[i]>>16&3], (int)(cns[i]>>16&0xffff)>>2);
	free(seqs);
	printf("//\n");
	fflush(stdout);
	g_vpos_shift += vpos;
}

static void update_vpos(int vpos, nseq_t *hash)
{
	khint_t k;
	for (k = 0; k < kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			rseq_t *p = &kh_val(hash, k);
			if (p->vpos < vpos) kh_del(64, hash, k); // TODO: if rseq_t::seq is allocated dynamically, free it
			else p->vpos -= vpos;
		}
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
		int i, j, c, cnt[4], tmp, dophase = 1;
		if (tid < 0) break;
		if (tid != lasttid) { // change of chromosome
			g_vpos_shift = 0;
			if (lasttid >= 0) phase(h->target_name[lasttid], vpos, cns, seqs);
			lasttid = tid;
			vpos = 0;
		}
		cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
		// check if there is a variant
		for (i = 0; i < n; ++i) {
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
		// add the variant
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
				dophase = 0;
			} else { // absent
				r->beg = p->b->core.pos;
				r->end = bam_calend(&p->b->core, bam1_cigar(p->b));
				r->vpos = vpos, r->vlen = 1, r->seq[0] = c;
			}
		}
		if (dophase) {
			phase(h->target_name[tid], vpos, cns, seqs);
			update_vpos(vpos, seqs);
			cns[0] = cns[vpos];
			vpos = 0;
		}
		++vpos;
	}
	if (tid >= 0) phase(h->target_name[tid], vpos, cns, seqs);
	bam_header_destroy(h);
	bam_plp_destroy(iter);
	bam_close(fp);
	kh_destroy(64, seqs);
	free(cns);
	return 0;
}
