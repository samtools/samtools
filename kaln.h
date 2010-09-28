/* The MIT License

   Copyright (c) 2003-2006, 2008, 2009 by Heng Li <lh3@live.co.uk>

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

#ifndef LH3_KALN_H_
#define LH3_KALN_H_

#include <stdint.h>

#define MINOR_INF -1073741823

typedef struct {
	int gap_open;
	int gap_ext;
	int gap_end_open;
	int gap_end_ext;

	int *matrix;
	int row;
	int band_width;
} ka_param_t;

typedef struct {
	float d, e;
	int bw;
} ka_probpar_t;

#ifdef __cplusplus
extern "C" {
#endif

	uint32_t *ka_global_core(uint8_t *seq1, int len1, uint8_t *seq2, int len2, const ka_param_t *ap,
							 int *_score, int *n_cigar);
	int ka_prob_glocal(const uint8_t *_ref, int l_ref, const uint8_t *_query, int l_query, const uint8_t *iqual,
					   const ka_probpar_t *c, int *state, uint8_t *q);

#ifdef __cplusplus
}
#endif

extern ka_param_t ka_param_blast; /* = { 5, 2, 5, 2, aln_sm_blast, 5, 50 }; */
extern ka_probpar_t ka_probpar_def; /* { 0.0001, 0.1, 10 } */

#endif
