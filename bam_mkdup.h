/* The MIT License
 *
 * Copyright (c) 2014 Genome Research Limited.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Contact: Martin Pollard <mp15@sanger.ac.uk> */

#ifndef BAM_MKDUP_H
#define BAM_MKDUP_H

struct possig;
typedef struct possig possig_t;

KHASH_DECLARE(sig, possig_t, char)

typedef struct state {
	samFile* fin;
	samFile* fout;
	bam_hdr_t* hin;
	bam_hdr_t* hout;
	int so;
	kh_sig_t* hash; // hash of seen positions could replace with bloom filter?
} state_t;

bool process_namesorted(const state_t* state);
bool process_coordsorted(/* HACK:const*/ state_t* state, const char* BIG_DIRTY_HACK);

static inline bool is_unprocessable( bam1_t* read )
{
	// if read is unpaired, unmapped, has unmapped mate, is secondary, qcfail or supplimentary
	return ((read->core.flag&(BAM_FUNMAP|BAM_FMUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FSUPPLEMENTARY)) != 0
			|| (read->core.flag&BAM_FPAIRED) != BAM_FPAIRED
			|| read->core.tid == -1
			|| read->core.pos == -1);
}

#endif // !defined(BAM_MKDUP_H)