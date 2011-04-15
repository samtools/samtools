/* The MIT License

   Copyright (c) 2008, 2010 by Attractive Chaos <attractor@live.co.uk>

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

/* Hooke-Jeeves algorithm for nonlinear minimization
 
   Based on the pseudocodes by Bell and Pike (CACM 9(9):684-685), and
   the revision by Tomlin and Smith (CACM 12(11):637-638). Both of the
   papers are comments on Kaupe's Algorithm 178 "Direct Search" (ACM
   6(6):313-314). The original algorithm was designed by Hooke and
   Jeeves (ACM 8:212-229). This program is further revised according to
   Johnson's implementation at Netlib (opt/hooke.c).
 
   Hooke-Jeeves algorithm is very simple and it works quite well on a
   few examples. However, it might fail to converge due to its heuristic
   nature. A possible improvement, as is suggested by Johnson, may be to
   choose a small r at the beginning to quickly approach to the minimum
   and a large r at later step to hit the minimum.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kmin.h"

static double __kmin_hj_aux(kmin_f func, int n, double *x1, void *data, double fx1, double *dx, int *n_calls)
{
	int k, j = *n_calls;
	double ftmp;
	for (k = 0; k != n; ++k) {
		x1[k] += dx[k];
		ftmp = func(n, x1, data); ++j;
		if (ftmp < fx1) fx1 = ftmp;
		else { /* search the opposite direction */
			dx[k] = 0.0 - dx[k];
			x1[k] += dx[k] + dx[k];
			ftmp = func(n, x1, data); ++j;
			if (ftmp < fx1) fx1 = ftmp;
			else x1[k] -= dx[k]; /* back to the original x[k] */
		}
	}
	*n_calls = j;
	return fx1; /* here: fx1=f(n,x1) */
}

double kmin_hj(kmin_f func, int n, double *x, void *data, double r, double eps, int max_calls)
{
	double fx, fx1, *x1, *dx, radius;
	int k, n_calls = 0;
	x1 = (double*)calloc(n, sizeof(double));
	dx = (double*)calloc(n, sizeof(double));
	for (k = 0; k != n; ++k) { /* initial directions, based on MGJ */
		dx[k] = fabs(x[k]) * r;
		if (dx[k] == 0) dx[k] = r;
	}
	radius = r;
	fx1 = fx = func(n, x, data); ++n_calls;
	for (;;) {
		memcpy(x1, x, n * sizeof(double)); /* x1 = x */
		fx1 = __kmin_hj_aux(func, n, x1, data, fx, dx, &n_calls);
		while (fx1 < fx) {
			for (k = 0; k != n; ++k) {
				double t = x[k];
				dx[k] = x1[k] > x[k]? fabs(dx[k]) : 0.0 - fabs(dx[k]);
				x[k] = x1[k];
				x1[k] = x1[k] + x1[k] - t;
			}
			fx = fx1;
			if (n_calls >= max_calls) break;
			fx1 = func(n, x1, data); ++n_calls;
			fx1 = __kmin_hj_aux(func, n, x1, data, fx1, dx, &n_calls);
			if (fx1 >= fx) break;
			for (k = 0; k != n; ++k)
				if (fabs(x1[k] - x[k]) > .5 * fabs(dx[k])) break;
			if (k == n) break;
		}
		if (radius >= eps) {
			if (n_calls >= max_calls) break;
			radius *= r;
			for (k = 0; k != n; ++k) dx[k] *= r;
		} else break; /* converge */
	}
	free(x1); free(dx);
	return fx1;
}

// I copied this function somewhere several years ago with some of my modifications, but I forgot the source.
double kmin_brent(kmin1_f func, double a, double b, void *data, double tol, double *xmin)
{
	double bound, u, r, q, fu, tmp, fa, fb, fc, c;
	const double gold1 = 1.6180339887;
	const double gold2 = 0.3819660113;
	const double tiny = 1e-20;
	const int max_iter = 100;

	double e, d, w, v, mid, tol1, tol2, p, eold, fv, fw;
	int iter;

	fa = func(a, data); fb = func(b, data);
	if (fb > fa) { // swap, such that f(a) > f(b)
		tmp = a; a = b; b = tmp;
		tmp = fa; fa = fb; fb = tmp;
	}
	c = b + gold1 * (b - a), fc = func(c, data); // golden section extrapolation
	while (fb > fc) {
		bound = b + 100.0 * (c - b); // the farthest point where we want to go
		r = (b - a) * (fb - fc);
		q = (b - c) * (fb - fa);
		if (fabs(q - r) < tiny) { // avoid 0 denominator
			tmp = q > r? tiny : 0.0 - tiny;
		} else tmp = q - r;
		u = b - ((b - c) * q - (b - a) * r) / (2.0 * tmp); // u is the parabolic extrapolation point
		if ((b > u && u > c) || (b < u && u < c)) { // u lies between b and c
			fu = func(u, data);
			if (fu < fc) { // (b,u,c) bracket the minimum
				a = b; b = u; fa = fb; fb = fu;
				break;
			} else if (fu > fb) { // (a,b,u) bracket the minimum
				c = u; fc = fu;
				break;
			}
			u = c + gold1 * (c - b); fu = func(u, data); // golden section extrapolation
		} else if ((c > u && u > bound) || (c < u && u < bound)) { // u lies between c and bound
			fu = func(u, data);
			if (fu < fc) { // fb > fc > fu
				b = c; c = u; u = c + gold1 * (c - b);
				fb = fc; fc = fu; fu = func(u, data);
			} else { // (b,c,u) bracket the minimum
				a = b; b = c; c = u;
				fa = fb; fb = fc; fc = fu;
				break;
			}
		} else if ((u > bound && bound > c) || (u < bound && bound < c)) { // u goes beyond the bound
			u = bound; fu = func(u, data);
		} else { // u goes the other way around, use golden section extrapolation
			u = c + gold1 * (c - b); fu = func(u, data);
		}
		a = b; b = c; c = u;
		fa = fb; fb = fc; fc = fu;
	}
	if (a > c) u = a, a = c, c = u; // swap

	// now, a<b<c, fa>fb and fb<fc, move on to Brent's algorithm
	e = d = 0.0;
	w = v = b; fv = fw = fb;
	for (iter = 0; iter != max_iter; ++iter) {
		mid = 0.5 * (a + c);
		tol2 = 2.0 * (tol1 = tol * fabs(b) + tiny);
		if (fabs(b - mid) <= (tol2 - 0.5 * (c - a))) {
			*xmin = b; return fb; // found
		}
		if (fabs(e) > tol1) {
			// related to parabolic interpolation
			r = (b - w) * (fb - fv);
			q = (b - v) * (fb - fw);
			p = (b - v) * q - (b - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) p = 0.0 - p;
			else q = 0.0 - q;
			eold = e; e = d;
			if (fabs(p) >= fabs(0.5 * q * eold) || p <= q * (a - b) || p >= q * (c - b)) {
				d = gold2 * (e = (b >= mid ? a - b : c - b));
			} else {
				d = p / q; u = b + d; // actual parabolic interpolation happens here
				if (u - a < tol2 || c - u < tol2)
					d = (mid > b)? tol1 : 0.0 - tol1;
			}
		} else d = gold2 * (e = (b >= mid ? a - b : c - b)); // golden section interpolation
		u = fabs(d) >= tol1 ? b + d : b + (d > 0.0? tol1 : -tol1);
		fu = func(u, data);
		if (fu <= fb) { // u is the minimum point so far
			if (u >= b) a = b;
			else c = b;
			v = w; w = b; b = u; fv = fw; fw = fb; fb = fu;
		} else { // adjust (a,c) and (u,v,w)
			if (u < b) a = u;
			else c = u;
			if (fu <= fw || w == b) {
				v = w; w = u;
				fv = fw; fw = fu;
			} else if (fu <= fv || v == b || v == w) {
				v = u; fv = fu;
			}
		}
	}
	*xmin = b;
	return fb;
}
