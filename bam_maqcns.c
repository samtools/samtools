#include <math.h>
#include <assert.h>
#include "bam.h"
#include "bam_maqcns.h"
#include "ksort.h"
#include "errmod.h"
#include "kaln.h"
KSORT_INIT_GENERIC(uint32_t)

#define INDEL_WINDOW_SIZE 50
#define INDEL_EXT_DEP 0.9

typedef struct __bmc_aux_t {
	int max;
	uint32_t *info;
	uint16_t *info16;
	errmod_t *em;
} bmc_aux_t;

typedef struct {
	float esum[4], fsum[4];
	uint32_t c[4];
} glf_call_aux_t;

/*
  P(<b1,b2>) = \theta \sum_{i=1}^{N-1} 1/i
  P(D|<b1,b2>) = \sum_{k=1}^{N-1} p_k 1/2 [(k/N)^n_2(1-k/N)^n_1 + (k/N)^n1(1-k/N)^n_2]
  p_k = 1/k / \sum_{i=1}^{N-1} 1/i
 */
static void cal_het(bam_maqcns_t *aa)
{
	int k, n1, n2;
	double sum_harmo; // harmonic sum
	double poly_rate;

	free(aa->lhet);
	aa->lhet = (double*)calloc(256 * 256, sizeof(double));
	sum_harmo = 0.0;
	for (k = 1; k <= aa->n_hap - 1; ++k)
		sum_harmo += 1.0 / k;
	for (n1 = 0; n1 < 256; ++n1) {
		for (n2 = 0; n2 < 256; ++n2) {
			long double sum = 0.0;
			double lC = aa->errmod == BAM_ERRMOD_SOAP? 0 : lgamma(n1+n2+1) - lgamma(n1+1) - lgamma(n2+1);
			for (k = 1; k <= aa->n_hap - 1; ++k) {
				double pk = 1.0 / k / sum_harmo;
				double log1 = log((double)k/aa->n_hap);
				double log2 = log(1.0 - (double)k/aa->n_hap);
				sum += pk * 0.5 * (expl(log1*n2) * expl(log2*n1) + expl(log1*n1) * expl(log2*n2));
			}
			aa->lhet[n1<<8|n2] = lC + logl(sum);
		}
	}
	poly_rate = aa->het_rate * sum_harmo;
	aa->q_r = -4.343 * log(2.0 * poly_rate / (1.0 - poly_rate));
}

/** initialize the helper structure */
static void cal_coef(bam_maqcns_t *aa)
{
	int k, n, q;
	long double sum_a[257], b[256], q_c[256], tmp[256], fk2[256];
	double *lC;

	if (aa->errmod == BAM_ERRMOD_MAQ2) return; // no need to do the following
	// aa->lhet will be allocated and initialized 
	free(aa->fk); free(aa->coef);
	aa->coef = 0;
	aa->fk = (double*)calloc(256, sizeof(double));
	aa->fk[0] = fk2[0] = 1.0;
	for (n = 1; n != 256; ++n) {
		aa->fk[n] = pow(aa->theta, n) * (1.0 - aa->eta) + aa->eta;
		fk2[n] = aa->fk[n>>1]; // this is an approximation, assuming reads equally likely come from both strands
	}
	if (aa->errmod == BAM_ERRMOD_SOAP) return;
	aa->coef = (double*)calloc(256*256*64, sizeof(double));
	lC = (double*)calloc(256 * 256, sizeof(double));
	for (n = 1; n != 256; ++n)
		for (k = 1; k <= n; ++k)
			lC[n<<8|k] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
	for (q = 1; q != 64; ++q) {
		double e = pow(10.0, -q/10.0);
		double le = log(e);
		double le1 = log(1.0-e);
		for (n = 1; n != 256; ++n) {
			double *coef = aa->coef + (q<<16|n<<8);
			sum_a[n+1] = 0.0;
			for (k = n; k >= 0; --k) { // a_k = \sum_{i=k}^n C^n_k \epsilon^k (1-\epsilon)^{n-k}
				sum_a[k] = sum_a[k+1] + expl(lC[n<<8|k] + k*le + (n-k)*le1);
				b[k] = sum_a[k+1] / sum_a[k];
				if (b[k] > 0.99) b[k] = 0.99;
			}
			for (k = 0; k != n; ++k) // log(\bar\beta_{nk}(\bar\epsilon)^{f_k})
				q_c[k] = -4.343 * fk2[k] * logl(b[k] / e);
			for (k = 1; k != n; ++k) q_c[k] += q_c[k-1]; // \prod_{i=0}^k c_i
			for (k = 0; k <= n; ++k) { // powl() in 64-bit mode seems broken on my Mac OS X 10.4.9
				tmp[k] = -4.343 * logl(1.0 - expl(fk2[k] * logl(b[k])));
				coef[k] = (k? q_c[k-1] : 0) + tmp[k]; // this is the final c_{nk}
			}
		}
	}
	free(lC);
}

bam_maqcns_t *bam_maqcns_init()
{
	bam_maqcns_t *bm;
	bm = (bam_maqcns_t*)calloc(1, sizeof(bam_maqcns_t));
	bm->aux = (bmc_aux_t*)calloc(1, sizeof(bmc_aux_t));
	bm->het_rate = 0.001;
	bm->theta = 0.83f;
	bm->n_hap = 2;
	bm->eta = 0.03;
	bm->cap_mapQ = 60;
	bm->min_baseQ = 13;
	return bm;
}

void bam_maqcns_prepare(bam_maqcns_t *bm)
{
	if (bm->errmod == BAM_ERRMOD_MAQ2) bm->aux->em = errmod_init(1. - bm->theta);
	cal_coef(bm); cal_het(bm);
}

void bam_maqcns_destroy(bam_maqcns_t *bm)
{
	if (bm == 0) return;
	free(bm->lhet); free(bm->fk); free(bm->coef); free(bm->aux->info); free(bm->aux->info16);
	if (bm->aux->em) errmod_destroy(bm->aux->em);
	free(bm->aux); free(bm);
}

glf1_t *bam_maqcns_glfgen(int _n, const bam_pileup1_t *pl, uint8_t ref_base, bam_maqcns_t *bm)
{
	glf_call_aux_t *b = 0;
	int i, j, k, w[8], c, n;
	glf1_t *g = (glf1_t*)calloc(1, sizeof(glf1_t));
	float p[16], min_p = 1e30;
	uint64_t rms;

	g->ref_base = ref_base;
	if (_n == 0) return g;

	// construct aux array
	if (bm->aux->max < _n) {
		bm->aux->max = _n;
		kroundup32(bm->aux->max);
		bm->aux->info = (uint32_t*)realloc(bm->aux->info, 4 * bm->aux->max);
		bm->aux->info16 = (uint16_t*)realloc(bm->aux->info16, 2 * bm->aux->max);
	}
	for (i = n = 0, rms = 0; i < _n; ++i) {
		const bam_pileup1_t *p = pl + i;
		uint32_t q, x = 0, qq;
		uint16_t y = 0;
		if (p->is_del || p->is_refskip || (p->b->core.flag&BAM_FUNMAP)) continue;
		q = (uint32_t)bam1_qual(p->b)[p->qpos];
		if (q < bm->min_baseQ) continue;
		x |= (uint32_t)bam1_strand(p->b) << 18 | q << 8 | p->b->core.qual;
		y |= bam1_strand(p->b)<<4;
		if (p->b->core.qual < q) q = p->b->core.qual;
		c = p->b->core.qual < bm->cap_mapQ? p->b->core.qual : bm->cap_mapQ;
		rms += c * c;
		x |= q << 24;
		y |= q << 5;
		qq = bam1_seqi(bam1_seq(p->b), p->qpos);
		q = bam_nt16_nt4_table[qq? qq : ref_base];
		if (!p->is_del && !p->is_refskip && q < 4) x |= 1 << 21 | q << 16, y |= q;
		bm->aux->info16[n] = y;
		bm->aux->info[n++] = x;
	}
	rms = (uint8_t)(sqrt((double)rms / n) + .499);
	if (bm->errmod == BAM_ERRMOD_MAQ2) {
		errmod_cal(bm->aux->em, n, 4, bm->aux->info16, p);
		goto goto_glf;
	}
	ks_introsort(uint32_t, n, bm->aux->info);
	// generate esum and fsum
	b = (glf_call_aux_t*)calloc(1, sizeof(glf_call_aux_t));
	for (k = 0; k != 8; ++k) w[k] = 0;
	for (j = n - 1; j >= 0; --j) { // calculate esum and fsum
		uint32_t info = bm->aux->info[j];
		if (info>>24 < 4 && (info>>8&0x3f) != 0) info = 4<<24 | (info&0xffffff);
		k = info>>16&7;
		if (info>>24 > 0) {
			b->esum[k&3] += bm->fk[w[k]] * (info>>24);
			b->fsum[k&3] += bm->fk[w[k]];
			if (w[k] < 0xff) ++w[k];
			++b->c[k&3];
		}
	}
	// rescale ->c[]
	for (j = c = 0; j != 4; ++j) c += b->c[j];
	if (c > 255) {
		for (j = 0; j != 4; ++j) b->c[j] = (int)(254.0 * b->c[j] / c + 0.5);
		for (j = c = 0; j != 4; ++j) c += b->c[j];
	}
	if (bm->errmod == BAM_ERRMOD_MAQ) {
		// generate likelihood
		for (j = 0; j != 4; ++j) {
			// homozygous
			float tmp1, tmp3;
			int tmp2, bar_e;
			for (k = 0, tmp1 = tmp3 = 0.0, tmp2 = 0; k != 4; ++k) {
				if (j == k) continue;
				tmp1 += b->esum[k]; tmp2 += b->c[k]; tmp3 += b->fsum[k];
			}
			if (tmp2) {
				bar_e = (int)(tmp1 / tmp3 + 0.5);
				if (bar_e < 4) bar_e = 4; // should not happen
				if (bar_e > 63) bar_e = 63;
				p[j<<2|j] = tmp1 + bm->coef[bar_e<<16|c<<8|tmp2];
			} else p[j<<2|j] = 0.0; // all the bases are j
			// heterozygous
			for (k = j + 1; k < 4; ++k) {
				for (i = 0, tmp2 = 0, tmp1 = tmp3 = 0.0; i != 4; ++i) {
					if (i == j || i == k) continue;
					tmp1 += b->esum[i]; tmp2 += b->c[i]; tmp3 += b->fsum[i];
				}
				if (tmp2) {
					bar_e = (int)(tmp1 / tmp3 + 0.5);
					if (bar_e < 4) bar_e = 4;
					if (bar_e > 63) bar_e = 63;
					p[j<<2|k] = p[k<<2|j] = -4.343 * bm->lhet[b->c[j]<<8|b->c[k]] + tmp1 + bm->coef[bar_e<<16|c<<8|tmp2];
				} else p[j<<2|k] = p[k<<2|j] = -4.343 * bm->lhet[b->c[j]<<8|b->c[k]]; // all the bases are either j or k
			}
			//
			for (k = 0; k != 4; ++k)
				if (p[j<<2|k] < 0.0) p[j<<2|k] = 0.0;
		}

		{ // fix p[k<<2|k]
			float max1, max2, min1, min2;
			int max_k, min_k;
			max_k = min_k = -1;
			max1 = max2 = -1.0; min1 = min2 = 1e30;
			for (k = 0; k < 4; ++k) {
				if (b->esum[k] > max1) {
					max2 = max1; max1 = b->esum[k]; max_k = k;
				} else if (b->esum[k] > max2) max2 = b->esum[k];
			}
			for (k = 0; k < 4; ++k) {
				if (p[k<<2|k] < min1) {
					min2 = min1; min1 = p[k<<2|k]; min_k = k;
				} else if (p[k<<2|k] < min2) min2 = p[k<<2|k];
			}
			if (max1 > max2 && (min_k != max_k || min1 + 1.0 > min2))
				p[max_k<<2|max_k] = min1 > 1.0? min1 - 1.0 : 0.0;
		}
	} else if (bm->errmod == BAM_ERRMOD_SOAP) { // apply the SOAP model
		// generate likelihood
		for (j = 0; j != 4; ++j) {
			float tmp;
			// homozygous
			for (k = 0, tmp = 0.0; k != 4; ++k)
				if (j != k) tmp += b->esum[k];
			p[j<<2|j] = tmp;
			// heterozygous
			for (k = j + 1; k < 4; ++k) {
				for (i = 0, tmp = 0.0; i != 4; ++i)
					if (i != j && i != k) tmp += b->esum[i];
				p[j<<2|k] = p[k<<2|j] = -4.343 * bm->lhet[b->c[j]<<8|b->c[k]] + tmp;
			}
		}
	}

goto_glf:
	// convert necessary information to glf1_t
	g->ref_base = ref_base; g->max_mapQ = rms;
	g->depth = n > 16777215? 16777215 : n;
	for (j = 0; j != 4; ++j)
		for (k = j; k < 4; ++k)
			if (p[j<<2|k] < min_p) min_p = p[j<<2|k];
	g->min_lk = min_p > 255.0? 255 : (int)(min_p + 0.5);
	for (j = c = 0; j != 4; ++j)
		for (k = j; k < 4; ++k)
			g->lk[c++] = p[j<<2|k]-min_p > 255.0? 255 : (int)(p[j<<2|k]-min_p + 0.5);

	free(b);
	return g;
}

uint32_t glf2cns(const glf1_t *g, int q_r)
{
	int i, j, k, p[10], ref4;
	uint32_t x = 0;
	ref4 = bam_nt16_nt4_table[g->ref_base];
	for (i = k = 0; i < 4; ++i)
		for (j = i; j < 4; ++j) {
			int prior = (i == ref4 && j == ref4? 0 : i == ref4 || j == ref4? q_r : q_r + 3);
			p[k] = (g->lk[k] + prior)<<4 | i<<2 | j;
			++k;
		}
	for (i = 1; i < 10; ++i) // insertion sort
		for (j = i; j > 0 && p[j] < p[j-1]; --j)
			k = p[j], p[j] = p[j-1], p[j-1] = k;
	x = (1u<<(p[0]&3) | 1u<<(p[0]>>2&3)) << 28; // the best genotype
	x |= (uint32_t)g->max_mapQ << 16; // rms mapQ
	x |= ((p[1]>>4) - (p[0]>>4) < 256? (p[1]>>4) - (p[0]>>4) : 255) << 8; // consensus Q
	for (k = 0; k < 10; ++k)
		if ((p[k]&0xf) == (ref4<<2|ref4)) break;
	if (k == 10) k = 9;
	x |= (p[k]>>4) - (p[0]>>4) < 256? (p[k]>>4) - (p[0]>>4) : 255; // snp Q
	return x;
}

uint32_t bam_maqcns_call(int n, const bam_pileup1_t *pl, bam_maqcns_t *bm)
{
	glf1_t *g;
	uint32_t x;
	if (n) {
		g = bam_maqcns_glfgen(n, pl, 0xf, bm);
		x = g->depth == 0? (0xfU<<28 | 0xfU<<24) : glf2cns(g, (int)(bm->q_r + 0.5));
		free(g);
	} else x = 0xfU<<28 | 0xfU<<24;
	return x;
}

/************** *****************/

bam_maqindel_opt_t *bam_maqindel_opt_init()
{
	bam_maqindel_opt_t *mi = (bam_maqindel_opt_t*)calloc(1, sizeof(bam_maqindel_opt_t));
	mi->q_indel = 40;
	mi->r_indel = 0.00015;
	mi->r_snp = 0.001;
	//
	mi->mm_penalty = 3;
	mi->indel_err = 4;
	mi->ambi_thres = 10;
	return mi;
}

void bam_maqindel_ret_destroy(bam_maqindel_ret_t *mir)
{
	if (mir == 0) return;
	free(mir->s[0]); free(mir->s[1]); free(mir);
}

int bam_tpos2qpos(const bam1_core_t *c, const uint32_t *cigar, int32_t tpos, int is_left, int32_t *_tpos)
{
	int k, x = c->pos, y = 0, last_y = 0;
	*_tpos = c->pos;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		int l = cigar[k] >> BAM_CIGAR_SHIFT;
		if (op == BAM_CMATCH) {
			if (c->pos > tpos) return y;
			if (x + l > tpos) {
				*_tpos = tpos;
				return y + (tpos - x);
			}
			x += l; y += l;
			last_y = y;
		} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
			if (x + l > tpos) {
				*_tpos = is_left? x : x + l;
				return y;
			}
			x += l;
		}
	}
	*_tpos = x;
	return last_y;
}

#define MINUS_CONST 0x10000000

bam_maqindel_ret_t *bam_maqindel(int n, int pos, const bam_maqindel_opt_t *mi, const bam_pileup1_t *pl, const char *ref,
								 int _n_types, int *_types)
{
	int i, j, n_types, *types, left, right, max_rd_len = 0;
	bam_maqindel_ret_t *ret = 0;
	// if there is no proposed indel, check if there is an indel from the alignment
	if (_n_types == 0) {
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pl + i;
			if (!(p->b->core.flag&BAM_FUNMAP) && p->indel != 0) break;
		}
		if (i == n) return 0; // no indel
	}
	{ // calculate how many types of indels are available (set n_types and types)
		int m;
		uint32_t *aux;
		aux = (uint32_t*)calloc(n + _n_types + 1, 4);
		m = 0;
		aux[m++] = MINUS_CONST; // zero indel is always a type
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pl + i;
			if (!(p->b->core.flag&BAM_FUNMAP) && p->indel != 0)
				aux[m++] = MINUS_CONST + p->indel;
			j = bam_cigar2qlen(&p->b->core, bam1_cigar(p->b));
			if (j > max_rd_len) max_rd_len = j;
		}
		if (_n_types) // then also add this to aux[]
			for (i = 0; i < _n_types; ++i)
				if (_types[i]) aux[m++] = MINUS_CONST + _types[i];
		ks_introsort(uint32_t, m, aux);
		// squeeze out identical types
		for (i = 1, n_types = 1; i < m; ++i)
			if (aux[i] != aux[i-1]) ++n_types;
		types = (int*)calloc(n_types, sizeof(int));
		j = 0;
		types[j++] = aux[0] - MINUS_CONST; 
		for (i = 1; i < m; ++i) {
			if (aux[i] != aux[i-1])
				types[j++] = aux[i] - MINUS_CONST;
		}
		free(aux);
	}
	{ // calculate left and right boundary
		left = pos > INDEL_WINDOW_SIZE? pos - INDEL_WINDOW_SIZE : 0;
		right = pos + INDEL_WINDOW_SIZE;
		if (types[0] < 0) right -= types[0];
		// in case the alignments stand out the reference
		for (i = pos; i < right; ++i)
			if (ref[i] == 0) break;
		right = i;
	}
	{ // the core part
		char *ref2, *rs, *inscns = 0;
		int qr_snp, k, l, *score, *pscore, max_ins = types[n_types-1];
		qr_snp = (int)(-4.343 * log(mi->r_snp) + .499);
		if (max_ins > 0) { // get the consensus of inserted sequences
			int *inscns_aux = (int*)calloc(4 * n_types * max_ins, sizeof(int));
			// count occurrences
			for (i = 0; i < n_types; ++i) {
				if (types[i] <= 0) continue; // not insertion
				for (j = 0; j < n; ++j) {
					const bam_pileup1_t *p = pl + j;
					if (!(p->b->core.flag&BAM_FUNMAP) && p->indel == types[i]) {
						for (k = 1; k <= p->indel; ++k) {
							int c = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos + k)];
							if (c < 4) ++inscns_aux[i*max_ins*4 + (k-1)*4 + c];
						}
					}
				}
			}
			// construct the consensus of inserted sequence
			inscns = (char*)calloc(n_types * max_ins, sizeof(char));
			for (i = 0; i < n_types; ++i) {
				for (j = 0; j < types[i]; ++j) {
					int max = 0, max_k = -1, *ia = inscns_aux + i*max_ins*4 + j*4;
					for (k = 0; k < 4; ++k) {
						if (ia[k] > max) {
							max = ia[k];
							max_k = k;
						}
					}
					inscns[i*max_ins + j] = max? 1<<max_k : 15;
				}
			}
			free(inscns_aux);
		}
		// calculate score
		ref2 = (char*)calloc(right - left + types[n_types-1] + 2, 1);
		rs   = (char*)calloc(right - left + max_rd_len + types[n_types-1] + 2, 1);
		score = (int*)calloc(n_types * n, sizeof(int));
		pscore = (int*)calloc(n_types * n, sizeof(int));
		for (i = 0; i < n_types; ++i) {
			ka_param_t ap = ka_param_blast;
			ap.band_width = 2 * types[n_types - 1] + 2;
			ap.gap_end_ext = 0;
			// write ref2
			for (k = 0, j = left; j <= pos; ++j)
				ref2[k++] = bam_nt16_nt4_table[bam_nt16_table[(int)ref[j]]];
			if (types[i] <= 0) j += -types[i];
			else for (l = 0; l < types[i]; ++l)
					 ref2[k++] = bam_nt16_nt4_table[(int)inscns[i*max_ins + l]];
			if (types[0] < 0) { // mask deleted sequences
				int jj, tmp = types[i] >= 0? -types[0] : -types[0] + types[i];
				for (jj = 0; jj < tmp && j < right && ref[j]; ++jj, ++j)
					ref2[k++] = 4;
			}
			for (; j < right && ref[j]; ++j)
				ref2[k++] = bam_nt16_nt4_table[bam_nt16_table[(int)ref[j]]];
			if (j < right) right = j;
			// calculate score for each read
			for (j = 0; j < n; ++j) {
				const bam_pileup1_t *p = pl + j;
				int qbeg, qend, tbeg, tend;
				if (p->b->core.flag & BAM_FUNMAP) continue;
				qbeg = bam_tpos2qpos(&p->b->core, bam1_cigar(p->b), left,  0, &tbeg);
				qend = bam_tpos2qpos(&p->b->core, bam1_cigar(p->b), right, 1, &tend);
				assert(tbeg >= left);
				for (l = qbeg; l < qend; ++l)
					rs[l - qbeg] = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), l)];
				{
					int x, y, n_acigar, ps;
					uint32_t *acigar;
					ps = 0;
					if (tend - tbeg + types[i] <= 0) {
						score[i*n+j] = -(1<<20);
						pscore[i*n+j] = 1<<20;
						continue;
					}
					acigar = ka_global_core((uint8_t*)ref2 + tbeg - left, tend - tbeg + types[i], (uint8_t*)rs, qend - qbeg, &ap, &score[i*n+j], &n_acigar);
					x = tbeg - left; y = 0;
					for (l = 0; l < n_acigar; ++l) {
						int op = acigar[l]&0xf;
						int len = acigar[l]>>4;
						if (op == BAM_CMATCH) {
							int k;
							for (k = 0; k < len; ++k)
								if (ref2[x+k] != rs[y+k] && ref2[x+k] < 4)
									ps += bam1_qual(p->b)[y+k] < qr_snp? bam1_qual(p->b)[y+k] : qr_snp;
							x += len; y += len;
						} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
							if (op == BAM_CINS && l > 0 && l < n_acigar - 1) ps += mi->q_indel * len;
							y += len;
						} else if (op == BAM_CDEL) {
							if (l > 0 && l < n_acigar - 1) ps += mi->q_indel * len;
							x += len;
						}
					}
					pscore[i*n+j] = ps;
					/*if (1) { // for debugging only
						fprintf(stderr, "id=%d, pos=%d, type=%d, j=%d, score=%d, psore=%d, %d, %d, %d, %d, %d, ",
								j, pos+1, types[i], j, score[i*n+j], pscore[i*n+j], tbeg, tend, qbeg, qend, mi->q_indel);
						for (l = 0; l < n_acigar; ++l) fprintf(stderr, "%d%c", acigar[l]>>4, "MIDS"[acigar[l]&0xf]);
						fprintf(stderr, "\n");
						for (l = 0; l < tend - tbeg + types[i]; ++l) fputc("ACGTN"[ref2[l+tbeg-left]], stderr);
						fputc('\n', stderr);
						for (l = 0; l < qend - qbeg; ++l) fputc("ACGTN"[rs[l]], stderr);
						fputc('\n', stderr);
						}*/
					free(acigar);
				}
			}
		}
		{ // get final result
			int *sum, max1, max2, max1_i, max2_i;
			// pick up the best two score
			sum = (int*)calloc(n_types, sizeof(int));
			for (i = 0; i < n_types; ++i)
				for (j = 0; j < n; ++j)
					sum[i] += -pscore[i*n+j];
			max1 = max2 = -0x7fffffff; max1_i = max2_i = -1;
			for (i = 0; i < n_types; ++i) {
				if (sum[i] > max1) {
					max2 = max1; max2_i = max1_i; max1 = sum[i]; max1_i = i;
				} else if (sum[i] > max2) {
					max2 = sum[i]; max2_i = i;
				}
			}
			free(sum);
			// write ret
			ret = (bam_maqindel_ret_t*)calloc(1, sizeof(bam_maqindel_ret_t));
			ret->indel1 = types[max1_i]; ret->indel2 = types[max2_i];
			ret->s[0] = (char*)calloc(abs(ret->indel1) + 2, 1);
			ret->s[1] = (char*)calloc(abs(ret->indel2) + 2, 1);
			// write indel sequence
			if (ret->indel1 > 0) {
				ret->s[0][0] = '+';
				for (k = 0; k < ret->indel1; ++k)
					ret->s[0][k+1] = bam_nt16_rev_table[(int)inscns[max1_i*max_ins + k]];
			} else if (ret->indel1 < 0) {
				ret->s[0][0] = '-';
				for (k = 0; k < -ret->indel1 && ref[pos + k + 1]; ++k)
					ret->s[0][k+1] = ref[pos + k + 1];
			} else ret->s[0][0] = '*';
			if (ret->indel2 > 0) {
				ret->s[1][0] = '+';
				for (k = 0; k < ret->indel2; ++k)
					ret->s[1][k+1] = bam_nt16_rev_table[(int)inscns[max2_i*max_ins + k]];
			} else if (ret->indel2 < 0) {
				ret->s[1][0] = '-';
				for (k = 0; k < -ret->indel2 && ref[pos + k + 1]; ++k)
					ret->s[1][k+1] = ref[pos + k + 1];
			} else ret->s[1][0] = '*';
			// write count
			for (i = 0; i < n; ++i) {
				const bam_pileup1_t *p = pl + i;
				if (p->indel == ret->indel1) ++ret->cnt1;
				else if (p->indel == ret->indel2) ++ret->cnt2;
				else ++ret->cnt_anti;
			}
			{ // write gl[]
				int tmp, seq_err = 0;
				double x = 1.0;
				tmp = max1_i - max2_i;
				if (tmp < 0) tmp = -tmp;
				for (j = 0; j < tmp + 1; ++j) x *= INDEL_EXT_DEP;
				seq_err = mi->q_indel * (1.0 - x) / (1.0 - INDEL_EXT_DEP);
				ret->gl[0] = ret->gl[1] = 0;
				for (j = 0; j < n; ++j) {
					int s1 = pscore[max1_i*n + j], s2 = pscore[max2_i*n + j];
					//fprintf(stderr, "id=%d, %d, %d, %d, %d, %d\n", j, pl[j].b->core.pos+1, types[max1_i], types[max2_i], s1, s2);
					if (s1 > s2) ret->gl[0] += s1 - s2 < seq_err? s1 - s2 : seq_err;
					else ret->gl[1] += s2 - s1 < seq_err? s2 - s1 : seq_err;
				}
			}
			// write cnt_ref and cnt_ambi
			if (max1_i != 0 && max2_i != 0) {
				for (j = 0; j < n; ++j) {
					int diff1 = score[j] - score[max1_i * n + j];
					int diff2 = score[j] - score[max2_i * n + j];
					if (diff1 > 0 && diff2 > 0) ++ret->cnt_ref;
					else if (diff1 == 0 || diff2 == 0) ++ret->cnt_ambi;
				}
			}
		}
		free(score); free(pscore); free(ref2); free(rs); free(inscns);
	}
	{ // call genotype
		int q[3], qr_indel = (int)(-4.343 * log(mi->r_indel) + 0.5);
		int min1, min2, min1_i;
		q[0] = ret->gl[0] + (ret->s[0][0] != '*'? 0 : 0) * qr_indel;
		q[1] = ret->gl[1] + (ret->s[1][0] != '*'? 0 : 0) * qr_indel;
		q[2] = n * 3 + (ret->s[0][0] == '*' || ret->s[1][0] == '*'? 1 : 1) * qr_indel;
		min1 = min2 = 0x7fffffff; min1_i = -1;
		for (i = 0; i < 3; ++i) {
			if (q[i] < min1) {
				min2 = min1; min1 = q[i]; min1_i = i;
			} else if (q[i] < min2) min2 = q[i];
		}
		ret->gt = min1_i;
		ret->q_cns = min2 - min1;
		// set q_ref
		if (ret->gt < 2) ret->q_ref = (ret->s[ret->gt][0] == '*')? 0 : q[1-ret->gt] - q[ret->gt] - qr_indel - 3;
		else ret->q_ref = (ret->s[0][0] == '*')? q[0] - q[2] : q[1] - q[2];
		if (ret->q_ref < 0) ret->q_ref = 0;
	}
	free(types);
	return ret;
}
