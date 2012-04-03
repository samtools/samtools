/* The MIT License

   Copyright (c) 2003-2006, 2008-2010, by Heng Li <lh3lh3@live.co.uk>

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
#include "kprobaln.h"

/*****************************************
 * Probabilistic banded glocal alignment *
 *****************************************/

#define EI .25
#define EM .33333333333

static float g_qual2prob[256];

#define set_u(u, b, i, k) { int x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }

kpa_par_t kpa_par_def = { 0.001, 0.1, 10 };
kpa_par_t kpa_par_alt = { 0.0001, 0.01, 10 };

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
int kpa_glocal(const uint8_t *_ref, int l_ref, const uint8_t *_query, int l_query, const uint8_t *iqual,
			   const kpa_par_t *c, int *state, uint8_t *q)
{
	double **f, **b = 0, *s, m[9], sI, sM, bI, bM, pb;
	float *qual, *_qual;
	const uint8_t *ref, *query;
	int bw, bw2, i, k, is_diff = 0, is_backward = 1, Pr;

    if ( l_ref<=0 || l_query<=0 ) return 0; // FIXME: this may not be an ideal fix, just prevents sefgault

	/*** initialization ***/
	is_backward = state && q? 1 : 0;
	ref = _ref - 1; query = _query - 1; // change to 1-based coordinate
	bw = l_ref > l_query? l_ref : l_query;
	if (bw > c->bw) bw = c->bw;
	if (bw < abs(l_ref - l_query)) bw = abs(l_ref - l_query);
	bw2 = bw * 2 + 1;
	// allocate the forward and backward matrices f[][] and b[][] and the scaling array s[]
	f = calloc(l_query+1, sizeof(void*));
	if (is_backward) b = calloc(l_query+1, sizeof(void*));
	for (i = 0; i <= l_query; ++i) {    // FIXME: this will lead in segfault for l_query==0
		f[i] = calloc(bw2 * 3 + 6, sizeof(double)); // FIXME: this is over-allocated for very short seqs
		if (is_backward) b[i] = calloc(bw2 * 3 + 6, sizeof(double));
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
	bM = (1 - c->d) / l_ref; bI = c->d / l_ref; // (bM+bI)*l_ref==1
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
	{ // compute likelihood
		double p = 1., Pr1 = 0.;
		for (i = 0; i <= l_query + 1; ++i) {
			p *= s[i];
			if (p < 1e-100) Pr1 += -4.343 * log(p), p = 1.;
		}
		Pr1 += -4.343 * log(p * l_ref * l_query);
		Pr = (int)(Pr1 + .499);
		if (!is_backward) { // skip backward and MAP
			for (i = 0; i <= l_query; ++i) free(f[i]);
			free(f); free(s); free(_qual);
			return Pr;
		}
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
		fprintf(stderr, "(%.10lg,%.10lg) (%d,%d:%c,%c:%d) %lg\n", pb, sum, i-1, max_k>>2,
				"ACGT"[query[i]], "ACGT"[ref[(max_k>>2)+1]], max_k&3, max); // DEBUG
#endif
	}
	/*** free ***/
	for (i = 0; i <= l_query; ++i) {
		free(f[i]); free(b[i]);
	}
	free(f); free(b); free(s); free(_qual);
	return Pr;
}

#ifdef _MAIN
#include <unistd.h>
int main(int argc, char *argv[])
{
	uint8_t conv[256], *iqual, *ref, *query;
	int c, l_ref, l_query, i, q = 30, b = 10, P;
	while ((c = getopt(argc, argv, "b:q:")) >= 0) {
		switch (c) {
		case 'b': b = atoi(optarg); break;
		case 'q': q = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: %s [-q %d] [-b %d] <ref> <query>\n", argv[0], q, b); // example: acttc attc
		return 1;
	}
	memset(conv, 4, 256);
	conv['a'] = conv['A'] = 0; conv['c'] = conv['C'] = 1;
	conv['g'] = conv['G'] = 2; conv['t'] = conv['T'] = 3;
	ref = (uint8_t*)argv[optind]; query = (uint8_t*)argv[optind+1];
	l_ref = strlen((char*)ref); l_query = strlen((char*)query);
	for (i = 0; i < l_ref; ++i) ref[i] = conv[ref[i]];
	for (i = 0; i < l_query; ++i) query[i] = conv[query[i]];
	iqual = malloc(l_query);
	memset(iqual, q, l_query);
	kpa_par_def.bw = b;
	P = kpa_glocal(ref, l_ref, query, l_query, iqual, &kpa_par_alt, 0, 0);
	fprintf(stderr, "%d\n", P);
	free(iqual);
	return 0;
}
#endif
