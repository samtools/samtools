#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include <zlib.h>
#include "bam.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define MAX_VARS 256
#define FLIP_PENALTY 2
#define FLIP_THRES 4
#define MASK_THRES 3

#define FLAG_FIX_CHIMERA 0x1
#define FLAG_MASK_POOR   0x2
#define FLAG_LIST_EXCL   0x4

typedef struct {
	// configurations, initialized in the main function
	int flag, k, min_mapQ, min_varQ;
	// other global variables
	int vpos_shift;
	bamFile fp;
	char *pre;
	bamFile out[4];
	// alignment queue
	int n, m;
	bam1_t **b;
} phaseg_t;

typedef struct {
	int8_t seq[MAX_VARS]; // TODO: change to dynamic memory allocation!
	int vpos, beg, end;
	uint32_t vlen:16, single:1, flip:1, phase:1, phased:1;
} frag_t, *frag_p;

#define rseq_lt(a,b) ((a)->vpos < (b)->vpos)

#include "khash.h"
KHASH_SET_INIT_INT64(set64)
KHASH_MAP_INIT_INT64(64, frag_t)

typedef khash_t(64) nseq_t;

#include "ksort.h"
KSORT_INIT(rseq, frag_p, rseq_lt)

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

static int **count_all(int l, int vpos, nseq_t *hash)
{
	khint_t k;
	int i, j, **cnt;
	uint8_t *seq;
	seq = calloc(l, 1);
	cnt = calloc(vpos, sizeof(void*));
	for (i = 0; i < vpos; ++i) cnt[i] = calloc(1<<l, sizeof(int));
	for (k = 0; k < kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			frag_t *p = &kh_val(hash, k);
			if (p->vpos >= vpos || p->single) continue; // out of region; or singleton
			if (p->vlen == 1) { // such reads should be flagged as deleted previously if everything is right
				p->single = 1;
				continue;
			}
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

// phasing
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

// phase each fragment
static uint64_t *fragphase(int vpos, const int8_t *path, nseq_t *hash, int flip)
{
	khint_t k;
	uint64_t *pcnt;
	uint32_t *left, *rght, max;
	left = rght = 0; max = 0;
	pcnt = calloc(vpos, 8);
	for (k = 0; k < kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			int i, c[2];
			frag_t *f = &kh_val(hash, k);
			if (f->vpos >= vpos) continue;
			// get the phase
			c[0] = c[1] = 0;
			for (i = 0; i < f->vlen; ++i) {
				if (f->seq[i] == 0) continue;
				++c[f->seq[i] == path[f->vpos + i] + 1? 0 : 1];
			}
			f->phase = c[0] > c[1]? 0 : 1;
			f->phased = c[0] == c[1]? 0 : 1;
			// fix chimera
			f->flip = 0;
			if (flip && c[0] >= 3 && c[1] >= 3) {
				int sum[2], m, mi, md;
				if (f->vlen > max) { // enlarge the array
					max = f->vlen;
					kroundup32(max);
					left = realloc(left, max * 4);
					rght = realloc(rght, max * 4);
				}
				for (i = 0, sum[0] = sum[1] = 0; i < f->vlen; ++i) { // get left counts
					if (f->seq[i]) {
						int c = f->phase? 2 - f->seq[i] : f->seq[i] - 1;
						++sum[c == path[f->vpos + i]? 0 : 1];
					}
					left[i] = sum[1]<<16 | sum[0];
				}
				for (i = f->vlen - 1, sum[0] = sum[1] = 0; i >= 0; --i) { // get right counts
					if (f->seq[i]) {
						int c = f->phase? 2 - f->seq[i] : f->seq[i] - 1;
						++sum[c == path[f->vpos + i]? 0 : 1];
					}
					rght[i] = sum[1]<<16 | sum[0];
				}
				// find the best flip point
				for (i = m = 0, mi = -1, md = -1; i < f->vlen - 1; ++i) {
					int a[2];
					a[0] = (left[i]&0xffff) + (rght[i+1]>>16&0xffff) - (rght[i+1]&0xffff) * FLIP_PENALTY;
					a[1] = (left[i]>>16&0xffff) + (rght[i+1]&0xffff) - (rght[i+1]>>16&0xffff) * FLIP_PENALTY;
					if (a[0] > a[1]) {
						if (a[0] > m) m = a[0], md = 0, mi = i;
					} else {
						if (a[1] > m) m = a[1], md = 1, mi = i;
					}
				}
				if (m - c[0] >= FLIP_THRES && m - c[1] >= FLIP_THRES) { // then flip
					f->flip = 1;
					if (md == 0) { // flip the tail
						for (i = mi + 1; i < f->vlen; ++i)
							if (f->seq[i] == 1) f->seq[i] = 2;
							else if (f->seq[i] == 2) f->seq[i] = 1;
					} else { // flip the head
						for (i = 0; i <= mi; ++i)
							if (f->seq[i] == 1) f->seq[i] = 2;
							else if (f->seq[i] == 2) f->seq[i] = 1;
					}
				}
			}
			// update pcnt[]
			if (!f->single) {
				for (i = 0; i < f->vlen; ++i) {
					int c;
					if (f->seq[i] == 0) continue;
					c = f->phase? 2 - f->seq[i] : f->seq[i] - 1;
					if (c == path[f->vpos + i]) {
						if (f->phase == 0) ++pcnt[f->vpos + i];
						else pcnt[f->vpos + i] += 1ull<<32;
					} else {
						if (f->phase == 0) pcnt[f->vpos + i] += 1<<16;
						else pcnt[f->vpos + i] += 1ull<<48;
					}
				}
			}
		}
	}
	free(left); free(rght);
	return pcnt;
}

static uint64_t *genmask(int vpos, const uint64_t *pcnt, int *_n)
{
	int i, max = 0, score = 0, max_i = -1, m = 0, n = 0, beg = 0;
	uint64_t *list = 0;
	for (i = 0; i < vpos; ++i) {
		uint64_t x = pcnt[i];
		int c = (x>>16&0xffff) + (x>>48&0xffff);
		int pre = score;
		//printf("%d\t%d\t%d\t%d\t%d\n", i, c, score, max, beg);
		score += c - 1;
		if ((x>>16&0xffff) > (x&0xffff)) score += (x>>16&0xffff) - (x&0xffff); // further penalty when there are more out-of-phase alleles
		if ((x>>48&0xffff) > (x>>32&0xffff)) score += (x>>48&0xffff) - (x>>32&0xffff);
		if (score < 0) score = 0;
		if (pre == 0 && score > 0) beg = i; // change from zero to non-zero
		if ((i == vpos - 1 || score == 0) && max >= MASK_THRES) {
			if (n == m) {
				m = m? m<<1 : 4;
				list = realloc(list, m * 8);
			}
			list[n++] = (uint64_t)beg<<32 | max_i;
			i = max_i; // reset i to max_i
			score = 0;
		} else if (score > max) max = score, max_i = i;
		if (score == 0) max = 0;
	}
	*_n = n;
	return list;
}

// trim heading and tailing ambiguous bases; mark deleted and remove sequence
static int clean_seqs(int vpos, nseq_t *hash)
{
	khint_t k;
	int ret = 0;
	for (k = 0; k < kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			frag_t *f = &kh_val(hash, k);
			int beg, end, i;
			if (f->vpos >= vpos) {
				ret = 1;
				continue;
			}
			for (i = 0; i < f->vlen; ++i)
				if (f->seq[i] != 0) break;
			beg = i;
			for (i = f->vlen - 1; i >= 0; --i)
				if (f->seq[i] != 0) break;
			end = i + 1;
			if (end - beg <= 0) kh_del(64, hash, k);
			else {
				if (beg != 0) memmove(f->seq, f->seq + beg, end - beg);
				f->vpos += beg; f->vlen = end - beg;
				f->single = f->vlen == 1? 1 : 0;
			}
		}
	}
	return ret;
}

// drop variants in specified regions; update cns and hash at the same time
static int dropreg(int vpos, int n_masked, const uint64_t *mask, const int8_t *path, uint64_t *cns, nseq_t *hash)
{
	int8_t *flt;
	int i, k, *map;
	khint_t j;
	flt = calloc(vpos, 1);
	map = calloc(vpos, sizeof(int));
	// get the list of sites to be filtered
	for (i = 0; i < n_masked; ++i)
		for (k = mask[i]>>32; k <= (int)mask[i]; ++k)
			flt[k] = 1;
	// generate map[]
	for (i = k = 0; i < vpos; ++i) {
		if (flt[i]) map[i] = -1;
		else map[i] = k++;
	}
	// filter hash
	for (j = 0; j < kh_end(hash); ++j) {
		if (kh_exist(hash, j)) {
			frag_t *s = &kh_val(hash, j);
			int new_vpos = -1;
			if (s->vpos >= vpos || s->single) continue;
			for (i = k = 0; i < s->vlen; ++i) {
				if (new_vpos < 0 && flt[s->vpos + i] == 0) new_vpos = s->vpos + i;
				if (flt[s->vpos + i] == 0)
					s->seq[k++] = s->seq[i];
			}
			if (k == 0) kh_del(64, hash, j); // no SNP
			else s->vlen = k, s->vpos = map[new_vpos], s->single = k==1? 1 : 0;
		}
	}
	// filter cns
	for (i = k = 0; i < vpos; ++i)
		if (flt[i] == 0) cns[k++] = cns[i];
	free(flt); free(map);
	return k;
}

static void dump_aln(phaseg_t *g, int min_pos, const nseq_t *hash)
{
	int i, is_flip;
	is_flip = (drand48() < 0.5);
	for (i = 0; i < g->n; ++i) {
		int end, which;
		uint64_t key;
		khint_t k;
		bam1_t *b = g->b[i];
		key = X31_hash_string(bam1_qname(b));
		end = bam_calend(&b->core, bam1_cigar(b));
		if (end > min_pos) break;
		k = kh_get(64, hash, key);
		if (k == kh_end(hash)) which = 3;
		else {
			frag_t *f = &kh_val(hash, k);
			if (f->phased && f->flip) which = 2;
			else if (f->phased == 0) which = 3;
			else which = f->phase;
			if (which < 2 && is_flip) which = 1 - which; // increase the randomness
		}
		bam_write1(g->out[which], b);
		bam_destroy1(b);
		g->b[i] = 0;
	}
	memmove(g->b, g->b + i, (g->n - i) * sizeof(void*));
	g->n -= i;
}

static int phase(phaseg_t *g, const char *chr, int vpos, uint64_t *cns, nseq_t *hash)
{
	int i, j, n_seqs = kh_size(hash), ori_vpos = vpos, n_masked = 0, min_pos;
	khint_t k;
	frag_t **seqs;
	int8_t *path;
	uint64_t *pcnt = 0, *regmask = 0;

	if (vpos < 2) return 0;
	i = clean_seqs(vpos, hash); // i is true if hash has an element with its vpos >= vpos
	min_pos = i? cns[vpos]>>32 : 0x7fffffff;
	{ // phase
		int **cnt;
		cnt = count_all(g->k, vpos, hash);
		path = dynaprog(g->k, vpos, cnt);
		for (i = 0; i < vpos; ++i) free(cnt[i]);
		free(cnt);
		if (g->flag & FLAG_MASK_POOR) {
			uint64_t *mask = 0;
			pcnt = fragphase(vpos, path, hash, 0); // do not fix chimeras during masking
			mask = genmask(vpos, pcnt, &n_masked);
			free(pcnt);
			regmask = calloc(n_masked, 8);
			for (i = 0; i < n_masked; ++i)
				regmask[i] = cns[mask[i]>>32]>>32<<32 | cns[(uint32_t)mask[i]]>>32;
			if ((vpos = dropreg(vpos, n_masked, mask, path, cns, hash)) < ori_vpos) {
				int last_i;
				free(path);
				clean_seqs(vpos, hash);
				cnt = count_all(g->k, vpos, hash);
				path = dynaprog(g->k, vpos, cnt);
				for (i = 1, last_i = 0; i < vpos; ++i) {
					int *c = cnt[i], j;
					for (j = 0; j < 1<<g->k && c[j] == 0; ++j);
					if (j == 1<<g->k) {
						printf("BL\t%s\t%d\t%d\n", chr, (int)(cns[last_i]>>32) + 1, (int)(cns[i-1]>>32) + 1);
						last_i = i;
					}
				}
				printf("BL\t%s\t%d\t%d\n", chr, (int)(cns[last_i]>>32) + 1, (int)(cns[vpos-1]>>32) + 1);
				for (i = 0; i < vpos; ++i) free(cnt[i]);
				free(cnt);
			} else printf("BL\t%s\t%d\t%d\n", chr, (int)(cns[0]>>32) + 1, (int)(cns[vpos-1]>>32) + 1);
			free(mask);
		}
		pcnt = fragphase(vpos, path, hash, g->flag & FLAG_FIX_CHIMERA);
	}
	if (regmask)
		for (i = 0; i < n_masked; ++i)
			printf("MK\t%s\t%d\t%d\n", chr, (int)(regmask[i]>>32) + 1, (int)regmask[i] + 1);
	for (i = 0; i < vpos; ++i) {
		uint64_t x = pcnt[i];
		int8_t c[2];
		c[0] = cns[i]&3; c[1] = cns[i]>>16&3;
		printf("VL\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%d\n", (int)(cns[i]>>32) + 1, i + g->vpos_shift + 1, "ACGT"[c[path[i]]], "ACGT"[c[1-path[i]]],
			(int)(x&0xffff), (int)(x>>16&0xffff), (int)(x>>32&0xffff), (int)(x>>48&0xffff));
	}
	free(path); free(pcnt); free(regmask);
	seqs = calloc(n_seqs, sizeof(void*));
	for (k = 0, i = 0; k < kh_end(hash); ++k) 
		if (kh_exist(hash, k) && kh_val(hash, k).vpos < vpos && !kh_val(hash, k).single)
			seqs[i++] = &kh_val(hash, k);
	n_seqs = i;
	ks_introsort_rseq(n_seqs, seqs);
	for (i = 0; i < n_seqs; ++i) {
		frag_t *s = seqs[i];
		printf("EV\t0\t%s\t%d\t40\t%dM\t*\t0\t0\t", chr, s->vpos + 1 + g->vpos_shift, s->vlen);
		for (j = 0; j < s->vlen; ++j) {
			uint32_t c = cns[s->vpos + j];
			if (s->seq[j] == 0) putchar('N');
			else putchar("ACGT"[s->seq[j] == 1? (c&3) : (c>>16&3)]);
		}
		printf("\t*\tYP:i:%d\tYF:i:%d\n", s->phase, s->flip);
	}
	free(seqs);
	printf("//\n");
	fflush(stdout);
	g->vpos_shift += vpos;
	dump_aln(g, min_pos, hash);
	return vpos;
}

static void update_vpos(int vpos, nseq_t *hash)
{
	khint_t k;
	for (k = 0; k < kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			frag_t *p = &kh_val(hash, k);
			if (p->vpos < vpos) kh_del(64, hash, k); // TODO: if frag_t::seq is allocated dynamically, free it
			else p->vpos -= vpos;
		}
	}
}

static nseq_t *shrink_hash(nseq_t *hash) // TODO: to implement
{
	return hash;
}

static int readaln(void *data, bam1_t *b)
{
	phaseg_t *g = (phaseg_t*)data;
	int ret;
	ret = bam_read1(g->fp, b);
	if (!(b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP)) && g->pre) {
		if (g->n == g->m) {
			g->m = g->m? g->m<<1 : 16;
			g->b = realloc(g->b, g->m * sizeof(void*));
		}
		g->b[g->n++] = bam_dup1(b);
	}
	return ret;
}

static khash_t(set64) *loadpos(const char *fn, bam_header_t *h)
{
	gzFile fp;
	kstream_t *ks;
	int ret, dret;
	kstring_t *str;
	khash_t(set64) *hash;

	hash = kh_init(set64);
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int tid = bam_get_tid(h, str->s);
		if (tid >= 0 && dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) >= 0) {
				uint64_t x = (uint64_t)tid<<32 | (atoi(str->s) - 1);
				kh_put(set64, hash, x, &ret);
			} else break;
		}
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (dret < 0) break;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return hash;
}

int main_phase(int argc, char *argv[])
{
	extern void bam_init_header_hash(bam_header_t *header);
	int c, tid, pos, vpos = 0, n, lasttid = -1, max_vpos = 0;
	const bam_pileup1_t *plp;
	bam_plp_t iter;
	bam_header_t *h;
	nseq_t *seqs;
	uint64_t *cns = 0;
	phaseg_t g;
	char *fn_list = 0;
	khash_t(set64) *set = 0;

	memset(&g, 0, sizeof(phaseg_t));
	g.flag = FLAG_FIX_CHIMERA | FLAG_MASK_POOR;
	g.min_varQ = 40; g.min_mapQ = 10; g.k = 11;
	while ((c = getopt(argc, argv, "eFMq:Q:k:b:l:")) >= 0) {
		switch (c) {
			case 'q': g.min_varQ = atoi(optarg); break;
			case 'Q': g.min_mapQ = atoi(optarg); break;
			case 'k': g.k = atoi(optarg); break;
			case 'F': g.flag &= ~FLAG_FIX_CHIMERA; break;
			case 'M': g.flag &= ~FLAG_MASK_POOR; break;
			case 'e': g.flag |= FLAG_LIST_EXCL; break;
			case 'b': g.pre = strdup(optarg); break;
			case 'l': fn_list = strdup(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools phase [options] <in.bam>\n\n");
		fprintf(stderr, "Options: -k INT    block length [%d]\n", g.k);
		fprintf(stderr, "         -b STR    prefix of BAMs to output [null]\n");
		fprintf(stderr, "         -q INT    min variant quality to call SNP [%d]\n", g.min_varQ);
		fprintf(stderr, "         -Q INT    min mapping quality [%d]\n", g.min_mapQ);
		fprintf(stderr, "         -l FILE   list of sites to phase [null]\n");
		fprintf(stderr, "         -F        do not attempt to fix chimeras\n");
		fprintf(stderr, "         -M        do not mask poorly phased regions\n");
		fprintf(stderr, "         -e        do not discover SNPs (effective with -l)\n");
		fprintf(stderr, "\n");
		return 1;
	}
	g.fp = bam_open(argv[optind], "r");
	h = bam_header_read(g.fp);
	if (fn_list) { // read the list of sites to phase
		bam_init_header_hash(h);
		set = loadpos(fn_list, h);
		free(fn_list);
	} else g.flag &= ~FLAG_LIST_EXCL;
	if (g.pre) { // open BAMs to write
		char *s = malloc(strlen(g.pre) + 20);
		strcpy(s, g.pre); strcat(s, ".0.bam"); g.out[0] = bam_open(s, "w");
		strcpy(s, g.pre); strcat(s, ".1.bam"); g.out[1] = bam_open(s, "w");
		strcpy(s, g.pre); strcat(s, ".chimera.bam"); g.out[2] = bam_open(s, "w");
		strcpy(s, g.pre); strcat(s, ".unphased.bam"); g.out[3] = bam_open(s, "w");
		for (c = 0; c <= 3; ++c) bam_header_write(g.out[c], h);
		free(s);
	}

	iter = bam_plp_init(readaln, &g);
	g.vpos_shift = 0;
	seqs = kh_init(64);
	while ((plp = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {
		int i, j, c, cnt[4], tmp, dophase = 1, in_set = 0;
		if (tid < 0) break;
		if (tid != lasttid) { // change of chromosome
			g.vpos_shift = 0;
			if (lasttid >= 0) {
				shrink_hash(seqs);
				phase(&g, h->target_name[lasttid], vpos, cns, seqs);
				update_vpos(0x7fffffff, seqs);
			}
			lasttid = tid;
			vpos = 0;
		}
		if (set && kh_get(set64, set, (uint64_t)tid<<32 | pos) != kh_end(set)) in_set = 1;
		cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
		// check if there is a variant
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = plp + i;
			uint8_t *seq = bam1_seq(p->b);
			uint8_t *qual = bam1_qual(p->b);
			if (p->is_del || p->is_refskip) continue;
			if (p->b->core.qual < g.min_mapQ) continue;
			cnt[(int)nt16_nt4_table[(int)bam1_seqi(seq, p->qpos)]] += qual[p->qpos];
		}
		for (i = 0; i < 4; ++i) {
			if (cnt[i] >= 1<<14) cnt[i] = (1<<14) - 1;
			cnt[i] = cnt[i]<<2 | i; // cnt[i] is 16-bit at most
		}
		for (i = 1; i < 4; ++i) // insertion sort
			for (j = i; j > 0 && cnt[j] > cnt[j-1]; --j)
				tmp = cnt[j], cnt[j] = cnt[j-1], cnt[j-1] = tmp;
		if (set && (g.flag&FLAG_LIST_EXCL) && !in_set) continue; // not in the list
		if (!in_set && cnt[1]>>2 <= g.min_varQ) continue; // not a variant
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
			frag_t *r;
			if (p->is_del || p->is_refskip) continue;
			if (p->b->core.qual < g.min_mapQ) continue;
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
				r->vpos = vpos, r->vlen = 1, r->seq[0] = c, r->single = 0;
			}
		}
		if (dophase) {
			seqs = shrink_hash(seqs);
			phase(&g, h->target_name[tid], vpos, cns, seqs);
			update_vpos(vpos, seqs);
			cns[0] = cns[vpos];
			vpos = 0;
		}
		++vpos;
	}
	if (tid >= 0) phase(&g, h->target_name[tid], vpos, cns, seqs);
	bam_header_destroy(h);
	bam_plp_destroy(iter);
	bam_close(g.fp);
	kh_destroy(64, seqs);
	kh_destroy(set64, set);
	free(cns);
	if (g.pre) {
		for (c = 0; c <= 3; ++c) bam_close(g.out[c]);
		free(g.pre); free(g.b);
	}
	return 0;
}
