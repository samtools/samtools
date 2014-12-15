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

#ifndef READ_VECTOR_H
#define READ_VECTOR_H

#include <htslib/sam.h>

#include <stdbool.h>
#include <inttypes.h>

struct read_vector;
typedef struct read_vector read_vector_t;

struct read_vector {
	bool orient:1;
	uint32_t tid:31;
	uint32_t pos;
};

static inline bool bam_to_read_vector(bam1_t* record, read_vector_t* part)
{
	if (record->core.tid > INT32_MAX) return false;
	part->orient = ((record->core.flag&BAM_FREVERSE) == BAM_FREVERSE);
	part->tid = record->core.tid;
    part->pos = record->core.pos;
	return true;
}

static inline bool read_vector_gt(read_vector_t left, read_vector_t right)
{
	if (left.tid > right.tid) {
		return true;
	} else if (right.tid > left.tid) {
		return false;
	} else {
		if (left.pos > right.pos) {
			return true;
		} else {
			return false;
		}
	}
}

#endif
