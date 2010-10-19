/* The MIT License

   Copyright (c) 2003-2006, 2008, 2009, by Heng Li <lh3lh3@gmail.com>

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "kaln.h"

#define FROM_M 0
#define FROM_I 1
#define FROM_D 2

typedef struct {
	int i, j;
	unsigned char ctype;
} path_t;

int aln_sm_blosum62[] = {
/*	 A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X */
	 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-4, 0,
	-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-4,-1,
	-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-4,-1,
	-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-4,-1,
	 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-4,-2,
	-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-4,-1,
	-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-4,-1,
	-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-4,-1,
	-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-4,-1,
	-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-1,
	-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-4,-1,
	-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-4,-1,
	-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-4,-2,
	 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,-4, 0,
	 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-4, 0,
	-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,
	-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-4,-1,
	 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-4,-1,
	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,
	 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-4,-1
};

int aln_sm_blast[] = {
	1, -3, -3, -3, -2,
	-3, 1, -3, -3, -2,
	-3, -3, 1, -3, -2,
	-3, -3, -3, 1, -2,
	-2, -2, -2, -2, -2
};

ka_param_t ka_param_blast = {  5,  2,   5, 2, aln_sm_blast, 5, 50 };
ka_param_t ka_param_aa2aa = { 10,  2,  10, 2, aln_sm_blosum62, 22, 50 };

static uint32_t *ka_path2cigar32(const path_t *path, int path_len, int *n_cigar)
{
	int i, n;
	uint32_t *cigar;
	unsigned char last_type;

	if (path_len == 0 || path == 0) {
		*n_cigar = 0;
		return 0;
	}

	last_type = path->ctype;
	for (i = n = 1; i < path_len; ++i) {
		if (last_type != path[i].ctype) ++n;
		last_type = path[i].ctype;
	}
	*n_cigar = n;
	cigar = (uint32_t*)calloc(*n_cigar, 4);

	cigar[0] = 1u << 4 | path[path_len-1].ctype;
	last_type = path[path_len-1].ctype;
	for (i = path_len - 2, n = 0; i >= 0; --i) {
		if (path[i].ctype == last_type) cigar[n] += 1u << 4;
		else {
			cigar[++n] = 1u << 4 | path[i].ctype;
			last_type = path[i].ctype;
		}
	}

	return cigar;
}

/***************************/
/* START OF common_align.c */
/***************************/

#define SET_INF(s) (s).M = (s).I = (s).D = MINOR_INF;

#define set_M(MM, cur, p, sc)							\
{														\
	if ((p)->M >= (p)->I) {								\
		if ((p)->M >= (p)->D) {							\
			(MM) = (p)->M + (sc); (cur)->Mt = FROM_M;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	} else {											\
		if ((p)->I > (p)->D) {							\
			(MM) = (p)->I + (sc); (cur)->Mt = FROM_I;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	}													\
}
#define set_I(II, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->I) {					\
		(cur)->It = FROM_M;								\
		(II) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->It = FROM_I;								\
		(II) = (p)->I - gap_ext;						\
	}													\
}
#define set_end_I(II, cur, p)							\
{														\
	if (gap_end_ext >= 0) {								\
		if ((p)->M - gap_end_open > (p)->I) {			\
			(cur)->It = FROM_M;							\
			(II) = (p)->M - gap_end_open - gap_end_ext;	\
		} else {										\
			(cur)->It = FROM_I;							\
			(II) = (p)->I - gap_end_ext;				\
		}												\
	} else set_I(II, cur, p);							\
}
#define set_D(DD, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->D) {					\
		(cur)->Dt = FROM_M;								\
		(DD) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->Dt = FROM_D;								\
		(DD) = (p)->D - gap_ext;						\
	}													\
}
#define set_end_D(DD, cur, p)							\
{														\
	if (gap_end_ext >= 0) {								\
		if ((p)->M - gap_end_open > (p)->D) {			\
			(cur)->Dt = FROM_M;							\
			(DD) = (p)->M - gap_end_open - gap_end_ext;	\
		} else {										\
			(cur)->Dt = FROM_D;							\
			(DD) = (p)->D - gap_end_ext;				\
		}												\
	} else set_D(DD, cur, p);							\
}

typedef struct {
	uint8_t Mt:3, It:2, Dt:2;
} dpcell_t;

typedef struct {
	int M, I, D;
} dpscore_t;

/***************************
 * banded global alignment *
 ***************************/
uint32_t *ka_global_core(uint8_t *seq1, int len1, uint8_t *seq2, int len2, const ka_param_t *ap, int *_score, int *n_cigar)
{
	int i, j;
	dpcell_t **dpcell, *q;
	dpscore_t *curr, *last, *s;
	int b1, b2, tmp_end;
	int *mat, end, max = 0;
	uint8_t type, ctype;
	uint32_t *cigar = 0;

	int gap_open, gap_ext, gap_end_open, gap_end_ext, b;
	int *score_matrix, N_MATRIX_ROW;

	/* initialize some align-related parameters. just for compatibility */
	gap_open = ap->gap_open;
	gap_ext = ap->gap_ext;
	gap_end_open = ap->gap_end_open;
	gap_end_ext = ap->gap_end_ext;
	b = ap->band_width;
	score_matrix = ap->matrix;
	N_MATRIX_ROW = ap->row;

	*n_cigar = 0;
	if (len1 == 0 || len2 == 0) return 0;

	/* calculate b1 and b2 */
	if (len1 > len2) {
		b1 = len1 - len2 + b;
		b2 = b;
	} else {
		b1 = b;
		b2 = len2 - len1 + b;
	}
	if (b1 > len1) b1 = len1;
	if (b2 > len2) b2 = len2;
	--seq1; --seq2;

	/* allocate memory */
	end = (b1 + b2 <= len1)? (b1 + b2 + 1) : (len1 + 1);
	dpcell = (dpcell_t**)malloc(sizeof(dpcell_t*) * (len2 + 1));
	for (j = 0; j <= len2; ++j)
		dpcell[j] = (dpcell_t*)malloc(sizeof(dpcell_t) * end);
	for (j = b2 + 1; j <= len2; ++j)
		dpcell[j] -= j - b2;
	curr = (dpscore_t*)malloc(sizeof(dpscore_t) * (len1 + 1));
	last = (dpscore_t*)malloc(sizeof(dpscore_t) * (len1 + 1));
	
	/* set first row */
	SET_INF(*curr); curr->M = 0;
	for (i = 1, s = curr + 1; i < b1; ++i, ++s) {
		SET_INF(*s);
		set_end_D(s->D, dpcell[0] + i, s - 1);
	}
	s = curr; curr = last; last = s;

	/* core dynamic programming, part 1 */
	tmp_end = (b2 < len2)? b2 : len2 - 1;
	for (j = 1; j <= tmp_end; ++j) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}
	/* last row for part 1, use set_end_D() instead of set_D() */
	if (j == len2 && b2 != len2 - 1) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_end_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
		++j;
	}

	/* core dynamic programming, part 2 */
	for (; j <= len2 - b2 + 1; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		end = j + b1 - 1;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}

	/* core dynamic programming, part 3 */
	for (; j < len2; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}
	/* last row */
	if (j == len2) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_end_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}

	*_score = last[len1].M;
	if (n_cigar) { /* backtrace */
		path_t *p, *path = (path_t*)malloc(sizeof(path_t) * (len1 + len2 + 2));
		i = len1; j = len2;
		q = dpcell[j] + i;
		s = last + len1;
		max = s->M; type = q->Mt; ctype = FROM_M;
		if (s->I > max) { max = s->I; type = q->It; ctype = FROM_I; }
		if (s->D > max) { max = s->D; type = q->Dt; ctype = FROM_D; }

		p = path;
		p->ctype = ctype; p->i = i; p->j = j; /* bug fixed 040408 */
		++p;
		do {
			switch (ctype) {
			case FROM_M: --i; --j; break;
			case FROM_I: --j; break;
			case FROM_D: --i; break;
			}
			q = dpcell[j] + i;
			ctype = type;
			switch (type) {
			case FROM_M: type = q->Mt; break;
			case FROM_I: type = q->It; break;
			case FROM_D: type = q->Dt; break;
			}
			p->ctype = ctype; p->i = i; p->j = j;
			++p;
		} while (i || j);
		cigar = ka_path2cigar32(path, p - path - 1, n_cigar);
		free(path);
	}

	/* free memory */
	for (j = b2 + 1; j <= len2; ++j)
		dpcell[j] += j - b2;
	for (j = 0; j <= len2; ++j)
		free(dpcell[j]);
	free(dpcell);
	free(curr); free(last);

	return cigar;
}

/*****************************************
 * Probabilistic banded glocal alignment *
 *****************************************/

static float g_qual2prob[256];

#define EI .25
#define EM .3333333333333
#define set_u(u, b, i, k) { int x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }

ka_probpar_t ka_probpar_def = { 0.001, 0.1, 10 };

/*
  The topology of the profile HMM:

           /\             /\        /\             /\
           I[1]           I[k-1]    I[k]           I[L]
            ^   \      \    ^    \   ^   \      \   ^
            |    \      \   |     \  |    \      \  |
    M[0]   M[1] -> ... -> M[k-1] -> M[k] -> ... -> M[L]   M[L+1]
                \      \/        \/      \/      /
                 \     /\        /\      /\     /
                       -> D[k-1] -> D[k] ->

   M[0] points to every {M,I}[k] and every {M,I}[k] points M[L+1].

   On input, _ref is the reference sequence and _query is the query
   sequence. Both are sequences of 0/1/2/3/4 where 4 stands for an
   ambiguous residue. iqual is the base quality. c sets the gap open
   probability, gap extension probability and band width.

   On output, state and q are arrays of length l_query. The higher 30
   bits give the reference position the query base is matched to and the
   lower two bits can be 0 (an alignment match) or 1 (an
   insertion). q[i] gives the phred scaled posterior probability of
   state[i] being wrong.
 */
int ka_prob_glocal(const uint8_t *_ref, int l_ref, const uint8_t *_query, int l_query, const uint8_t *iqual,
				   const ka_probpar_t *c, int *state, uint8_t *q)
{
	double **f, **b, *s, m[9], sI, sM, bI, bM, pb;
	float *qual, *_qual;
	const uint8_t *ref, *query;
	int bw, bw2, i, k, is_diff = 0;

	/*** initialization ***/
	ref = _ref - 1; query = _query - 1; // change to 1-based coordinate
	bw = l_ref > l_query? l_ref : l_query;
	if (bw > c->bw) bw = c->bw;
	if (bw < abs(l_ref - l_query)) bw = abs(l_ref - l_query);
	bw2 = bw * 2 + 1;
	// allocate the forward and backward matrices f[][] and b[][] and the scaling array s[]
	f = calloc(l_query+1, sizeof(void*));
	b = calloc(l_query+1, sizeof(void*));
	for (i = 0; i <= l_query; ++i) {
		f[i] = calloc(bw2 * 3 + 6, sizeof(double)); // FIXME: this is over-allocated for very short seqs
		b[i] = calloc(bw2 * 3 + 6, sizeof(double));
	}
	s = calloc(l_query+2, sizeof(double)); // s[] is the scaling factor to avoid underflow
	// initialize qual
	_qual = calloc(l_query, sizeof(float));
	if (g_qual2prob[0] == 0)
		for (i = 0; i < 256; ++i)
			g_qual2prob[i] = pow(10, -i/10.);
	for (i = 0; i < l_query; ++i) _qual[i] = g_qual2prob[iqual? iqual[i] : 30];
	qual = _qual - 1;
	// initialize transition probability
	sM = sI = 1. / (2 * l_query + 2); // the value here seems not to affect results; FIXME: need proof
	m[0*3+0] = (1 - c->d - c->d) * (1 - sM); m[0*3+1] = m[0*3+2] = c->d * (1 - sM);
	m[1*3+0] = (1 - c->e) * (1 - sI); m[1*3+1] = c->e * (1 - sI); m[1*3+2] = 0.;
	m[2*3+0] = 1 - c->e; m[2*3+1] = 0.; m[2*3+2] = c->e;
	bM = (1 - c->d) / l_query; bI = c->d / l_query; // (bM+bI)*l_query==1
	/*** forward ***/
	// f[0]
	set_u(k, bw, 0, 0);
	f[0][k] = s[0] = 1.;
	{ // f[1]
		double *fi = f[1], sum;
		int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1, _beg, _end;
		for (k = beg, sum = 0.; k <= end; ++k) {
			int u;
			double e = (ref[k] > 3 || query[1] > 3)? 1. : ref[k] == query[1]? 1. - qual[1] : qual[1] * EM;
			set_u(u, bw, 1, k);
			fi[u+0] = e * bM; fi[u+1] = EI * bI;
			sum += fi[u] + fi[u+1];
		}
		// rescale
		s[1] = sum;
		set_u(_beg, bw, 1, beg); set_u(_end, bw, 1, end); _end += 2;
		for (k = _beg; k <= _end; ++k) fi[k] /= sum;
	}
	// f[2..l_query]
	for (i = 2; i <= l_query; ++i) {
		double *fi = f[i], *fi1 = f[i-1], sum, qli = qual[i];
		int beg = 1, end = l_ref, x, _beg, _end;
		uint8_t qyi = query[i];
		x = i - bw; beg = beg > x? beg : x; // band start
		x = i + bw; end = end < x? end : x; // band end
		for (k = beg, sum = 0.; k <= end; ++k) {
			int u, v11, v01, v10;
			double e;
			e = (ref[k] > 3 || qyi > 3)? 1. : ref[k] == qyi? 1. - qli : qli * EM;
			set_u(u, bw, i, k); set_u(v11, bw, i-1, k-1); set_u(v10, bw, i-1, k); set_u(v01, bw, i, k-1);
			fi[u+0] = e * (m[0] * fi1[v11+0] + m[3] * fi1[v11+1] + m[6] * fi1[v11+2]);
			fi[u+1] = EI * (m[1] * fi1[v10+0] + m[4] * fi1[v10+1]);
			fi[u+2] = m[2] * fi[v01+0] + m[8] * fi[v01+2];
			sum += fi[u] + fi[u+1] + fi[u+2];
//			fprintf(stderr, "F (%d,%d;%d): %lg,%lg,%lg\n", i, k, u, fi[u], fi[u+1], fi[u+2]); // DEBUG
		}
		// rescale
		s[i] = sum;
		set_u(_beg, bw, i, beg); set_u(_end, bw, i, end); _end += 2;
		for (k = _beg, sum = 1./sum; k <= _end; ++k) fi[k] *= sum;
	}
	{ // f[l_query+1]
		double sum;
		for (k = 1, sum = 0.; k <= l_ref; ++k) {
			int u;
			set_u(u, bw, l_query, k);
			if (u < 3 || u >= bw2*3+3) continue;
		    sum += f[l_query][u+0] * sM + f[l_query][u+1] * sI;
		}
		s[l_query+1] = sum; // the last scaling factor
	}
	/*** backward ***/
	// b[l_query] (b[l_query+1][0]=1 and thus \tilde{b}[][]=1/s[l_query+1]; this is where s[l_query+1] comes from)
	for (k = 1; k <= l_ref; ++k) {
		int u;
		double *bi = b[l_query];
		set_u(u, bw, l_query, k);
		if (u < 3 || u >= bw2*3+3) continue;
		bi[u+0] = sM / s[l_query] / s[l_query+1]; bi[u+1] = sI / s[l_query] / s[l_query+1];
	}
	// b[l_query-1..1]
	for (i = l_query - 1; i >= 1; --i) {
		int beg = 1, end = l_ref, x, _beg, _end;
		double *bi = b[i], *bi1 = b[i+1], y = (i > 1), qli1 = qual[i+1];
		uint8_t qyi1 = query[i+1];
		x = i - bw; beg = beg > x? beg : x;
		x = i + bw; end = end < x? end : x;
		for (k = end; k >= beg; --k) {
			int u, v11, v01, v10;
			double e;
			set_u(u, bw, i, k); set_u(v11, bw, i+1, k+1); set_u(v10, bw, i+1, k); set_u(v01, bw, i, k+1);
			e = (k >= l_ref? 0 : (ref[k+1] > 3 || qyi1 > 3)? 1. : ref[k+1] == qyi1? 1. - qli1 : qli1 * EM) * bi1[v11];
			bi[u+0] = e * m[0] + EI * m[1] * bi1[v10+1] + m[2] * bi[v01+2]; // bi1[v11] has been foled into e.
			bi[u+1] = e * m[3] + EI * m[4] * bi1[v10+1];
			bi[u+2] = (e * m[6] + m[8] * bi[v01+2]) * y;
//			fprintf(stderr, "B (%d,%d;%d): %lg,%lg,%lg\n", i, k, u, bi[u], bi[u+1], bi[u+2]); // DEBUG
		}
		// rescale
		set_u(_beg, bw, i, beg); set_u(_end, bw, i, end); _end += 2;
		for (k = _beg, y = 1./s[i]; k <= _end; ++k) bi[k] *= y;
	}
	{ // b[0]
		int beg = 1, end = l_ref < bw + 1? l_ref : bw + 1;
		double sum = 0.;
		for (k = end; k >= beg; --k) {
			int u;
			double e = (ref[k] > 3 || query[1] > 3)? 1. : ref[k] == query[1]? 1. - qual[1] : qual[1] * EM;
			set_u(u, bw, 1, k);
			if (u < 3 || u >= bw2*3+3) continue;
		    sum += e * b[1][u+0] * bM + EI * b[1][u+1] * bI;
		}
		set_u(k, bw, 0, 0);
		pb = b[0][k] = sum / s[0]; // if everything works as is expected, pb == 1.0
	}
	is_diff = fabs(pb - 1.) > 1e-7? 1 : 0;
	/*** MAP ***/
	for (i = 1; i <= l_query; ++i) {
		double sum = 0., *fi = f[i], *bi = b[i], max = 0.;
		int beg = 1, end = l_ref, x, max_k = -1;
		x = i - bw; beg = beg > x? beg : x;
		x = i + bw; end = end < x? end : x;
		for (k = beg; k <= end; ++k) {
			int u;
			double z;
			set_u(u, bw, i, k);
			z = fi[u+0] * bi[u+0]; if (z > max) max = z, max_k = (k-1)<<2 | 0; sum += z;
			z = fi[u+1] * bi[u+1]; if (z > max) max = z, max_k = (k-1)<<2 | 1; sum += z;
		}
		max /= sum; sum *= s[i]; // if everything works as is expected, sum == 1.0
		if (state) state[i-1] = max_k;
		if (q) k = (int)(-4.343 * log(1. - max) + .499), q[i-1] = k > 100? 99 : k;
#ifdef _MAIN
		fprintf(stderr, "(%.10lg,%.10lg) (%d,%d:%d)~%lg\n", pb, sum, i-1, max_k>>2, max_k&3, max); // DEBUG
#endif
	}
	/*** free ***/
	for (i = 0; i <= l_query; ++i) {
		free(f[i]); free(b[i]);
	}
	free(f); free(b); free(s); free(_qual);
	return 0;
}

#ifdef _MAIN
int main()
{
	int l_ref = 5, l_query = 4;
	uint8_t *ref = (uint8_t*)"\0\1\3\3\1";
	uint8_t *query = (uint8_t*)"\0\3\3\1";
//	uint8_t *query = (uint8_t*)"\1\3\3\1";
	static uint8_t qual[4] = {20, 20, 20, 20};
	ka_prob_glocal(ref, l_ref, query, l_query, qual, &ka_probpar_def, 0, 0);
	return 0;
}
#endif
