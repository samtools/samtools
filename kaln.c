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

int aln_sm_qual[] = {
	  0, -23, -23, -23, 0,
	-23,   0, -23, -23, 0,
	-23, -23,   0, -23, 0,
	-23, -23, -23,   0, 0,
	  0,   0,   0,   0, 0
};

ka_param_t ka_param_blast = {  5,  2,   5, 2, aln_sm_blast, 5, 50 };
ka_param_t ka_param_aa2aa = { 10,  2,  10, 2, aln_sm_blosum62, 22, 50 };

ka_param2_t ka_param2_qual  = { 37, 11, 37, 11, 37, 11, 0, 0, aln_sm_qual, 5, 50 };

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
	uint8_t Mt:3, It:2, Dt:3;
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

	if (n_cigar) *n_cigar = 0;
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

typedef struct {
	int M, I, D;
} score_aux_t;

#define MINUS_INF -0x40000000

// matrix: len2 rows and len1 columns
int ka_global_score(const uint8_t *_seq1, int len1, const uint8_t *_seq2, int len2, const ka_param2_t *ap)
{
	
#define __score_aux(_p, _q0, _sc, _io, _ie, _do, _de) {					\
		int t1, t2;														\
		score_aux_t *_q;												\
		_q = _q0;														\
		_p->M = _q->M >= _q->I? _q->M : _q->I;							\
		_p->M = _p->M >= _q->D? _p->M : _q->D;							\
		_p->M += (_sc);													\
		++_q;      t1 = _q->M - _io - _ie; t2 = _q->I - _ie; _p->I = t1 >= t2? t1 : t2; \
		_q = _p-1; t1 = _q->M - _do - _de; t2 = _q->D - _de; _p->D = t1 >= t2? t1 : t2; \
	}

	int i, j, bw, scmat_size = ap->row, *scmat = ap->matrix, ret;
	const uint8_t *seq1, *seq2;
	score_aux_t *curr, *last, *swap;
	bw = abs(len1 - len2) + ap->band_width;
	i = len1 > len2? len1 : len2;
	if (bw > i + 1) bw = i + 1;
	seq1 = _seq1 - 1; seq2 = _seq2 - 1;
	curr = calloc(len1 + 2, sizeof(score_aux_t));
	last = calloc(len1 + 2, sizeof(score_aux_t));
	{ // the zero-th row
		int x, end = len1;
		score_aux_t *p;
		j = 0;
		x = j + bw; end = len1 < x? len1 : x; // band end
		p = curr;
		p->M = 0; p->I = p->D = MINUS_INF;
		for (i = 1, p = &curr[1]; i <= end; ++i, ++p)
			p->M = p->I = MINUS_INF, p->D = -(ap->edo + ap->ede * i);
		p->M = p->I = p->D = MINUS_INF;
		swap = curr; curr = last; last = swap;
	}
	for (j = 1; j < len2; ++j) {
		int x, beg = 0, end = len1, *scrow, col_end;
		score_aux_t *p;
		x = j - bw; beg =    0 > x?    0 : x; // band start
		x = j + bw; end = len1 < x? len1 : x; // band end
		if (beg == 0) { // from zero-th column
			p = curr;
			p->M = p->D = MINUS_INF; p->I = -(ap->eio + ap->eie * j);
			++beg; // then beg = 1
		}
		scrow = scmat + seq2[j] * scmat_size;
		if (end == len1) col_end = 1, --end;
		else col_end = 0;
		for (i = beg, p = &curr[beg]; i <= end; ++i, ++p)
			__score_aux(p, &last[i-1], scrow[(int)seq1[i]], ap->iio, ap->iie, ap->ido, ap->ide);
		if (col_end) {
			__score_aux(p, &last[i-1], scrow[(int)seq1[i]], ap->eio, ap->eie, ap->ido, ap->ide);
			++p;
		}
		p->M = p->I = p->D = MINUS_INF;
//		for (i = 0; i <= len1; ++i) printf("(%d,%d,%d) ", curr[i].M, curr[i].I, curr[i].D); putchar('\n');
		swap = curr; curr = last; last = swap;
	}
	{ // the last row
		int x, beg = 0, *scrow;
		score_aux_t *p;
		j = len2;
		x = j - bw; beg = 0 > x?    0 : x; // band start
		if (beg == 0) { // from zero-th column
			p = curr;
			p->M = p->D = MINUS_INF; p->I = -(ap->eio + ap->eie * j);
			++beg; // then beg = 1
		}
		scrow = scmat + seq2[j] * scmat_size;
		for (i = beg, p = &curr[beg]; i < len1; ++i, ++p)
			__score_aux(p, &last[i-1], scrow[(int)seq1[i]], ap->iio, ap->iie, ap->edo, ap->ede);
		__score_aux(p, &last[i-1], scrow[(int)seq1[i]], ap->eio, ap->eie, ap->edo, ap->ede);
//		for (i = 0; i <= len1; ++i) printf("(%d,%d,%d) ", curr[i].M, curr[i].I, curr[i].D); putchar('\n');
	}
	ret = curr[len1].M >= curr[len1].I? curr[len1].M : curr[len1].I;
	ret = ret >= curr[len1].D? ret : curr[len1].D;
	free(curr); free(last);
	return ret;
}

#ifdef _MAIN
int main(int argc, char *argv[])
{
//	int len1 = 35, len2 = 35;
//	uint8_t *seq1 = (uint8_t*)"\0\0\3\3\2\0\0\0\1\0\2\1\2\1\3\2\3\3\3\0\2\3\2\1\1\3\3\3\2\3\3\1\0\0\1";
//	uint8_t *seq2 = (uint8_t*)"\0\0\3\3\2\0\0\0\1\0\2\1\2\1\3\2\3\3\3\0\2\3\2\1\1\3\3\3\2\3\3\1\0\1\0";
	int len1 = 4, len2 = 4;
	uint8_t *seq1 = (uint8_t*)"\1\0\0\1";
	uint8_t *seq2 = (uint8_t*)"\1\0\1\0";
	int sc;
//	ka_global_core(seq1, 2, seq2, 1, &ka_param_qual, &sc, 0);
	sc = ka_global_score(seq1, len1, seq2, len2, &ka_param2_qual);
	printf("%d\n", sc);
	return 0;
}
#endif
