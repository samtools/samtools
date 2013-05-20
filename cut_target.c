#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "bam.h"
#include "errmod.h"
#include "faidx.h"

#define ERR_DEP 0.83

typedef struct {
	int e[2][3], p[2][2];
} score_param_t;

/* Note that although the two matrics have 10 parameters in total, only 4
 * (probably 3) are free.  Changing the scoring matrices in a sort of symmetric
 * way will not change the result. */
static score_param_t g_param = { {{0,0,0},{-4,1,6}}, {{0,-14000}, {0,0}} };

typedef struct {
	int min_baseQ, tid, max_bases;
	uint16_t *bases;
	bamFile fp;
	bam_header_t *h;
	char *ref;
	faidx_t *fai;
	errmod_t *em;
} ct_t;

static uint16_t gencns(ct_t *g, int n, const bam_pileup1_t *plp)
{
	int i, j, ret, tmp, k, sum[4], qual;
	float q[16];
	if (n > g->max_bases) { // enlarge g->bases
		g->max_bases = n;
		kroundup32(g->max_bases);
		g->bases = realloc(g->bases, g->max_bases * 2);
	}
	for (i = k = 0; i < n; ++i) {
		const bam_pileup1_t *p = plp + i;
		uint8_t *seq;
		int q, baseQ, b;
		if (p->is_refskip || p->is_del) continue;
		baseQ = bam1_qual(p->b)[p->qpos];
		if (baseQ < g->min_baseQ) continue;
		seq = bam1_seq(p->b);
		b = bam_nt16_nt4_table[bam1_seqi(seq, p->qpos)];
		if (b > 3) continue;
		q = baseQ < p->b->core.qual? baseQ : p->b->core.qual;
		if (q < 4) q = 4;
		if (q > 63) q = 63;
		g->bases[k++] = q<<5 | bam1_strand(p->b)<<4 | b;
	}
	if (k == 0) return 0;
	errmod_cal(g->em, k, 4, g->bases, q);
	for (i = 0; i < 4; ++i) sum[i] = (int)(q[i<<2|i] + .499) << 2 | i;
	for (i = 1; i < 4; ++i) // insertion sort
		for (j = i; j > 0 && sum[j] < sum[j-1]; --j)
			tmp = sum[j], sum[j] = sum[j-1], sum[j-1] = tmp;
	qual = (sum[1]>>2) - (sum[0]>>2);
	k = k < 256? k : 255;
	ret = (qual < 63? qual : 63) << 2 | (sum[0]&3);
	return ret<<8|k;
}

static void process_cns(bam_header_t *h, int tid, int l, uint16_t *cns)
{
	int i, f[2][2], *prev, *curr, *swap_tmp, s;
	uint8_t *b; // backtrack array
	b = calloc(l, 1);
	f[0][0] = f[0][1] = 0;
	prev = f[0]; curr = f[1];
	// fill the backtrack matrix
	for (i = 0; i < l; ++i) {
		int c = (cns[i] == 0)? 0 : (cns[i]>>8 == 0)? 1 : 2;
		int tmp0, tmp1;
		// compute f[0]
		tmp0 = prev[0] + g_param.e[0][c] + g_param.p[0][0]; // (s[i+1],s[i])=(0,0)
		tmp1 = prev[1] + g_param.e[0][c] + g_param.p[1][0]; // (0,1)
		if (tmp0 > tmp1) curr[0] = tmp0, b[i] = 0;
		else curr[0] = tmp1, b[i] = 1;
		// compute f[1]
		tmp0 = prev[0] + g_param.e[1][c] + g_param.p[0][1]; // (s[i+1],s[i])=(1,0)
		tmp1 = prev[1] + g_param.e[1][c] + g_param.p[1][1]; // (1,1)
		if (tmp0 > tmp1) curr[1] = tmp0, b[i] |= 0<<1;
		else curr[1] = tmp1, b[i] |= 1<<1;
		// swap
		swap_tmp = prev; prev = curr; curr = swap_tmp;
	}
	// backtrack
	s = prev[0] > prev[1]? 0 : 1;
	for (i = l - 1; i > 0; --i) {
		b[i] |= s<<2;
		s = b[i]>>s&1;
	}
	// print
	for (i = 0, s = -1; i <= l; ++i) {
		if (i == l || ((b[i]>>2&3) == 0 && s >= 0)) {
			if (s >= 0) {
				int j;
				printf("%s:%d-%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t", h->target_name[tid], s+1, i, h->target_name[tid], s+1, i-s);
				for (j = s; j < i; ++j) {
					int c = cns[j]>>8;
					if (c == 0) putchar('N');
					else putchar("ACGT"[c&3]);
				}
				putchar('\t');
				for (j = s; j < i; ++j)
					putchar(33 + (cns[j]>>8>>2));
				putchar('\n');
			}
			//if (s >= 0) printf("%s\t%d\t%d\t%d\n", h->target_name[tid], s, i, i - s);
			s = -1;
		} else if ((b[i]>>2&3) && s < 0) s = i;
	}
	free(b);
}

static int read_aln(void *data, bam1_t *b)
{
	extern int bam_prob_realn_core(bam1_t *b, const char *ref, int flag);
	ct_t *g = (ct_t*)data;
	int ret, len;
	ret = bam_read1(g->fp, b);
	if (ret >= 0 && g->fai && b->core.tid >= 0 && (b->core.flag&4) == 0) {
		if (b->core.tid != g->tid) { // then load the sequence
			free(g->ref);
			g->ref = fai_fetch(g->fai, g->h->target_name[b->core.tid], &len);
			g->tid = b->core.tid;
		}
		bam_prob_realn_core(b, g->ref, 1<<1|1);
	}
	return ret;
}

int main_cut_target(int argc, char *argv[])
{
	int c, tid, pos, n, lasttid = -1, lastpos = -1, l, max_l;
	const bam_pileup1_t *p;
	bam_plp_t plp;
	uint16_t *cns;
	ct_t g;

	memset(&g, 0, sizeof(ct_t));
	g.min_baseQ = 13; g.tid = -1;
	while ((c = getopt(argc, argv, "f:Q:i:o:0:1:2:")) >= 0) {
		switch (c) {
			case 'Q': g.min_baseQ = atoi(optarg); break; // quality cutoff
			case 'i': g_param.p[0][1] = -atoi(optarg); break; // 0->1 transition (in) PENALTY
			case '0': g_param.e[1][0] = atoi(optarg); break; // emission SCORE
			case '1': g_param.e[1][1] = atoi(optarg); break;
			case '2': g_param.e[1][2] = atoi(optarg); break;
			case 'f': g.fai = fai_load(optarg);
				if (g.fai == 0) fprintf(stderr, "[%s] fail to load the fasta index.\n", __func__);
				break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samtools targetcut [-Q minQ] [-i inPen] [-0 em0] [-1 em1] [-2 em2] [-f ref] <in.bam>\n");
		return 1;
	}
	l = max_l = 0; cns = 0;
	g.fp = strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	g.h = bam_header_read(g.fp);
	g.em = errmod_init(1. - ERR_DEP);
	plp = bam_plp_init(read_aln, &g);
	while ((p = bam_plp_auto(plp, &tid, &pos, &n)) != 0) {
		if (tid < 0) break;
		if (tid != lasttid) { // change of chromosome
			if (cns) process_cns(g.h, lasttid, l, cns);
			if (max_l < g.h->target_len[tid]) {
				max_l = g.h->target_len[tid];
				kroundup32(max_l);
				cns = realloc(cns, max_l * 2);
			}
			l = g.h->target_len[tid];
			memset(cns, 0, max_l * 2);
			lasttid = tid;
		}
		cns[pos] = gencns(&g, n, p);
		lastpos = pos;
	}
	process_cns(g.h, lasttid, l, cns);
	free(cns);
	bam_header_destroy(g.h);
	bam_plp_destroy(plp);
	bam_close(g.fp);
	if (g.fai) {
		fai_destroy(g.fai); free(g.ref);
	}
	errmod_destroy(g.em);
	free(g.bases);
	return 0;
}
