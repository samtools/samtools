#include <math.h>
#include "bam.h"
#include "bam_maqcns.h"
#include "ksort.h"
KSORT_INIT_GENERIC(uint32_t)

typedef struct __bmc_aux_t {
	int max;
	uint32_t *info;
} bmc_aux_t;

typedef struct {
	float esum[4], fsum[4];
	uint32_t c[4];
	uint32_t mapQ_max;
} glf_call_aux_t;

char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

/*
  P(<b1,b2>) = \theta \sum_{i=1}^{N-1} 1/i
  P(D|<b1,b2>) = \sum_{k=1}^{N-1} p_k 1/2 [(k/N)^n_2(1-k/N)^n_1 + (k/N)^n1(1-k/N)^n_2]
  p_k = i/k / \sum_{i=1}^{N-1} 1/i
 */
static void cal_het(bam_maqcns_t *aa)
{
	int k, n1, n2;
	double sum_harmo; // harmonic sum
	double poly_rate;
	double p1 = 0.0, p3 = 0.0; // just for testing

	free(aa->lhet);
	aa->lhet = (double*)calloc(256 * 256, sizeof(double));
	sum_harmo = 0.0;
	for (k = 1; k <= aa->n_hap - 1; ++k)
		sum_harmo += 1.0 / k;
	for (n1 = 0; n1 < 256; ++n1) {
		for (n2 = 0; n2 < 256; ++n2) {
			long double sum = 0.0;
			double lC = lgamma(n1+n2+1) - lgamma(n1+1) - lgamma(n2+1); // \binom{n1+n2}{n1}
			for (k = 1; k <= aa->n_hap - 1; ++k) {
				double pk = 1.0 / k / sum_harmo;
				double log1 = log((double)k/aa->n_hap);
				double log2 = log(1.0 - (double)k/aa->n_hap);
				sum += pk * 0.5 * (expl(log1*n2) * expl(log2*n1) + expl(log1*n1) * expl(log2*n2));
			}
			aa->lhet[n1<<8|n2] = lC + logl(sum);
			if (n1 == 17 && n2 == 3) p3 = lC + logl(expl(logl(0.5) * 20));
			if (n1 == 19 && n2 == 1) p1 = lC + logl(expl(logl(0.5) * 20));
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

	lC = (double*)calloc(256 * 256, sizeof(double));
	// aa->lhet will be allocated and initialized 
	free(aa->fk); free(aa->coef);
	aa->fk = (double*)calloc(256, sizeof(double));
	aa->coef = (double*)calloc(256*256*64, sizeof(double));
	aa->fk[0] = fk2[0] = 1.0;
	for (n = 1; n != 256; ++n) {
		aa->fk[n] = pow(aa->theta, n) * (1.0 - aa->eta) + aa->eta;
		fk2[n] = aa->fk[n>>1]; // this is an approximation, assuming reads equally likely come from both strands
	}
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
	bm->theta = 0.85;
	bm->n_hap = 2;
	bm->eta = 0.03;
	return bm;
}

void bam_maqcns_prepare(bam_maqcns_t *bm)
{
	cal_coef(bm); cal_het(bm);
}

void bam_maqcns_destroy(bam_maqcns_t *bm)
{
	if (bm == 0) return;
	free(bm->lhet); free(bm->fk); free(bm->coef); free(bm->aux->info);
	free(bm->aux); free(bm);
}

glf1_t *bam_maqcns_glfgen(int _n, const bam_pileup1_t *pl, uint8_t ref_base, bam_maqcns_t *bm)
{
	glf_call_aux_t *b;
	int i, j, k, w[8], c, n;
	glf1_t *g = (glf1_t*)calloc(1, sizeof(glf1_t));
	float p[16], min_p = 1e30;

	g->ref_base = ref_base;
	if (_n == 0) return g;

	// construct aux array
	if (bm->aux->max < _n) {
		bm->aux->max = _n;
		kroundup32(bm->aux->max);
		bm->aux->info = (uint32_t*)realloc(bm->aux->info, 4 * bm->aux->max);
	}
	for (i = n = 0; i < _n; ++i) {
		const bam_pileup1_t *p = pl + i;
		uint32_t q, x = 0;
		if (p->is_del || (p->b->core.flag&BAM_FUNMAP)) continue;
		q = (uint32_t)bam1_qual(p->b)[p->qpos];
		x |= (uint32_t)bam1_strand(p->b) << 18 | q << 8 | p->b->core.qual;
		if (p->b->core.qual < q) q = p->b->core.qual;
		x |= q << 24;
		q = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
		if (!p->is_del && q < 4) x |= 1 << 21 | q << 16;
		bm->aux->info[n++] = x;
	}
	ks_introsort(uint32_t, n, bm->aux->info);
	// generate esum and fsum
	b = (glf_call_aux_t*)calloc(1, sizeof(glf_call_aux_t));
	for (k = 0; k != 8; ++k) w[k] = 0;
	b->mapQ_max = 0;
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
		if (b->mapQ_max < (info&0x7f)) b->mapQ_max = info&0x7f;
	}
	// rescale ->c[]
	for (j = c = 0; j != 4; ++j) c += b->c[j];
	if (c > 255) {
		for (j = 0; j != 4; ++j) b->c[j] = (int)(254.0 * b->c[j] / c + 0.5);
		for (j = c = 0; j != 4; ++j) c += b->c[j];
	}
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

	// convert necessary information to glf1_t
	g->ref_base = ref_base; g->max_mapQ = b->mapQ_max;
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
	int i, j, k, tmp[16], min = 10000, min2 = 10000, min3 = 10000, min_g = -1, min_g2 = -1;
	uint32_t x = 0;
	for (i = k = 0; i < 4; ++i)
		for (j = i; j < 4; ++j) {
			tmp[j<<2|i] = -1;
			tmp[i<<2|j] = g->lk[k++] + (i == j? 0 : q_r);
		}
	for (i = 0; i < 16; ++i) {
		if (tmp[i] < 0) continue;
		if (tmp[i] < min) {
			min3 = min2; min2 = min; min = tmp[i]; min_g2 = min_g; min_g = i;
		} else if (tmp[i] < min2) {
			min3 = min2; min2 = tmp[i]; min_g2 = i;
		} else if (tmp[i] < min3) min3 = tmp[i];
	}
	x = min_g >= 0? (1U<<(min_g>>2&3) | 1U<<(min_g&3)) << 28 : 0xf << 28;
	x |= min_g2 >= 0? (1U<<(min_g2>>2&3) | 1U<<(min_g2&3)) << 24 : 0xf << 24;
	x |= (uint32_t)g->max_mapQ << 16;
	x |= min2 < 10000? (min2 - min < 256? min2 - min : 255) << 8 : 0xff << 8;
	x |= min2 < 10000 && min3 < 10000? (min3 - min2 < 256? min3 - min2 : 255) : 0xff;
	return x;
}

uint32_t bam_maqcns_call(int n, const bam_pileup1_t *pl, bam_maqcns_t *bm)
{
	glf1_t *g;
	uint32_t x;
	if (n) {
		g = bam_maqcns_glfgen(n, pl, 0xf, bm);
		x = glf2cns(g, (int)(bm->q_r + 0.5));
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

#define MINUS_CONST 0x10000000

bam_maqindel_ret_t *bam_maqindel(int n, int pos, const bam_maqindel_opt_t *mi, const bam_pileup1_t *pl, const char *ref)
{
	int i, j, n_types, *types, left, right;
	bam_maqindel_ret_t *ret = 0;
	for (i = 0; i < n; ++i) {
		const bam_pileup1_t *p = pl + i;
		if (!(p->b->core.flag&BAM_FUNMAP) && p->indel != 0) break;
	}
	if (i == n) return 0; // no indel
	{ // calculate how many types of indels are available (set n_types and types)
		int m;
		uint32_t *aux;
		aux = (uint32_t*)calloc(n+1, 4);
		m = 0;
		aux[m++] = MINUS_CONST; // zero indel is always a type
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pl + i;
			if (!(p->b->core.flag&BAM_FUNMAP) && p->indel != 0)
				aux[m++] = MINUS_CONST + p->indel;
		}
		ks_introsort(uint32_t, m, aux);
		n_types = 1;
		for (i = 1; i < m; ++i)
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
		bam_segreg_t seg;
		left = 0x7fffffff; right = 0;
		for (i = 0; i < n; ++i) {
			const bam_pileup1_t *p = pl + i;
			if (!(p->b->core.flag&BAM_FUNMAP)) {
				bam_segreg(pos, &p->b->core, bam1_cigar(p->b), &seg);
				if (seg.tbeg < left) left = seg.tbeg;
				if (seg.tend > right) right = seg.tend;
			}
		}
	}
	{ // the core part
		char *ref2, *inscns = 0;
		int k, l, *score, *pscore, max_ins = types[n_types-1];
		ref2 = (char*)calloc(right - left + types[n_types-1] + 2, 1);
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
		score = (int*)calloc(n_types * n, sizeof(int));
		pscore = (int*)calloc(n_types * n, sizeof(int));
		for (i = 0; i < n_types; ++i) {
			// write ref2
			for (k = 0, j = left; j <= pos; ++j)
				ref2[k++] = bam_nt16_table[(int)ref[j]];
			if (types[i] <= 0) j += -types[i];
			else for (l = 0; l < types[i]; ++l)
					 ref2[k++] = inscns[i*max_ins + l];
			for (; j < right && ref[j]; ++j)
				ref2[k++] = bam_nt16_table[(int)ref[j]];
			// calculate score for each read
			for (j = 0; j < n; ++j) {
				const bam_pileup1_t *p = pl + j;
				uint32_t *cigar;
				bam1_core_t *c = &p->b->core;
				int s, ps;
				bam_segreg_t seg;
				if (c->flag&BAM_FUNMAP) continue;
				cigar = bam1_cigar(p->b);
				bam_segreg(pos, c, cigar, &seg);
				for (ps = s = 0, l = seg.qbeg; c->pos + l < right && l < seg.qend; ++l) {
					int cq = bam1_seqi(bam1_seq(p->b), l), ct;
					ct = c->pos + l >= left? ref2[c->pos + l - left] : 15; // "<" should not happen if there is no bug
					if (cq < 15 && ct < 15) {
						s += cq == ct? 1 : -mi->mm_penalty;
						if (cq != ct) ps += bam1_qual(p->b)[l];
					}
				}
				score[i*n + j] = s; pscore[i*n + j] = ps;
				if (types[i] != 0) { // then try the other way to calculate the score
					for (ps = s = 0, l = seg.qbeg; c->pos + l + types[i] < right && l < seg.qend; ++l) {
						int cq = bam1_seqi(bam1_seq(p->b), l), ct;
						ct = c->pos + l + types[i] >= left? ref2[c->pos + l + types[i] - left] : 15;
						if (cq < 15 && ct < 15) {
							s += cq == ct? 1 : -mi->mm_penalty;
							if (cq != ct) ps += bam1_qual(p->b)[l];
						}
					}
				}
				if (score[i*n+j] < s) score[i*n+j] = s; // choose the higher of the two scores
				if (pscore[i*n+j] > ps) pscore[i*n+j] = ps;
				if (types[i] != 0) score[i*n+j] -= mi->indel_err;
				//printf("%d, %d, %d, %d\n", i, types[i], j, score[i*n+j]);
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
			for (j = 0; j < n; ++j) {
				if (score[max1_i*n+j] < 0 && score[max2_i*n+j] < 0) ++ret->cnt_anti;
				else {
					int diff = score[max1_i*n+j] - score[max2_i*n+j];
					if (diff > mi->ambi_thres) ++ret->cnt1;
					else if (diff < -mi->ambi_thres) ++ret->cnt2;
					else ++ret->cnt_ambi;
				}
			}
			// write gl[]
			ret->gl[0] = ret->gl[1] = 0;
			for (j = 0; j < n; ++j) {
				int s1 = pscore[max1_i*n + j], s2 = pscore[max2_i*n + j];
				ret->gl[0] += s1 < s2? 0 : s1 - s2 < mi->q_indel? s1 - s2 : mi->q_indel;
				ret->gl[1] += s2 < s1? 0 : s2 - s1 < mi->q_indel? s2 - s1 : mi->q_indel;
			}
		}
		free(score); free(pscore); free(ref2); free(inscns);
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
