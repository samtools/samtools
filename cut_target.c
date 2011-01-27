#include <unistd.h>
#include <stdlib.h>
#include "bam.h"

typedef struct {
	int e[2][3], p[2][2];
} score_param_t;

/* Note that although the two matrics have 10 parameters in total, there are
 * only 4 FREE parameters. Changing the scoring matrices in a sort of symmetric
 * way will not change the result. */
static score_param_t param = { {{0,0,0},{-4,1,6}}, {{0,-14000}, {0,0}} };

static int min_Q = 20;

static uint16_t gencns(int n, const bam_pileup1_t *plp)
{
	int i, j, ret, depth = 0, sum[4], tmp;
	sum[0] = sum[1] = sum[2] = sum[3] = 0;
	for (i = 0; i < n; ++i) {
		const bam_pileup1_t *p = plp + i;
		uint8_t *seq;
		int c, q;
		if (p->is_refskip) continue;
		++depth;
		if (p->is_del || p->b->core.qual < min_Q) continue;
		q = bam1_qual(p->b)[p->qpos];
		if (q < min_Q) continue;
		seq = bam1_seq(p->b);
		c = bam_nt16_nt4_table[bam1_seqi(seq, p->qpos)];
		if (c > 3) continue;
		sum[c] += q;
	}
	if (depth == 0) return 0;
	for (i = 0; i < 4; ++i) sum[i] = sum[i]<<2 | i;
	for (i = 1; i < 4; ++i) // insertion sort
		for (j = i; j > 0 && sum[j] > sum[j-1]; --j)
			tmp = sum[j], sum[j] = sum[j-1], sum[j-1] = tmp;
	depth = depth < 256? depth : 255;
	if (sum[0]>>2 == 0) return depth; 
	if (sum[1]>>2 > min_Q) ret = 61<<2 | (sum[0]&3); // a het
	else ret = (sum[0]>>2 < 60? sum[0]>>2 : 60) << 2 | (sum[0]&3);
	return depth<<8|ret;
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
		tmp0 = prev[0] + param.e[0][c] + param.p[0][0]; // (s[i+1],s[i])=(0,0)
		tmp1 = prev[1] + param.e[0][c] + param.p[1][0]; // (0,1)
		if (tmp0 > tmp1) curr[0] = tmp0, b[i] = 0;
		else curr[0] = tmp1, b[i] = 1;
		// compute f[1]
		tmp0 = prev[0] + param.e[1][c] + param.p[0][1]; // (s[i+1],s[i])=(1,0)
		tmp1 = prev[1] + param.e[1][c] + param.p[1][1]; // (1,1)
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
	// TODO: split overlapping or false joining fosmids

	// print
	for (i = 0, s = -1; i <= l; ++i) {
		if (i == l || ((b[i]>>2&3) == 0 && s >= 0)) {
			if (s >= 0) printf("%s\t%d\t%d\t%d\n", h->target_name[tid], s, i, i - s);
			s = -1;
		} else if ((b[i]>>2&3) && s < 0) s = i;
	}
	free(b);
}

int main_cut_target(int argc, char *argv[])
{
	bamFile fp;
	int c, tid, pos, n, lasttid = -1, lastpos = -1, l, max_l;
	const bam_pileup1_t *p;
	bam_plp_t plp;
	bam_header_t *h;
	uint16_t *cns;

	while ((c = getopt(argc, argv, "Q:i:o:0:1:2:")) >= 0) {
		switch (c) {
			case 'Q': min_Q = atoi(optarg); break; // quality cutoff
			case 'i': param.p[0][1] = -atoi(optarg); break; // 0->1 transition (in) PENALTY
			case '0': param.e[1][0] = atoi(optarg); break; // emission SCORE
			case '1': param.e[1][1] = atoi(optarg); break;
			case '2': param.e[1][2] = atoi(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samtools targetcut <in.bam>\n");
		return 1;
	}
	l = max_l = 0; cns = 0;
	fp = bam_open(argv[optind], "r");
	h = bam_header_read(fp);
	plp = bam_plp_init((bam_plp_auto_f)bam_read1, fp);
	while ((p = bam_plp_auto(plp, &tid, &pos, &n)) != 0) {
		if (tid < 0) break;
		if (tid != lasttid) { // change of chromosome
			if (cns) process_cns(h, lasttid, l, cns);
			if (max_l < h->target_len[tid]) {
				max_l = h->target_len[tid];
				kroundup32(max_l);
				cns = realloc(cns, max_l * 2);
			}
			l = h->target_len[tid];
			memset(cns, 0, max_l * 2);
			lasttid = tid;
		}
		cns[pos] = gencns(n, p);
		lastpos = pos;
	}
	process_cns(h, lasttid, l, cns);
	free(cns);
	bam_header_destroy(h);
	bam_plp_destroy(plp);
	bam_close(fp);
	return 0;
}
