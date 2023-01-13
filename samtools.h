/*  samtools.h -- utility routines.

    Copyright (C) 2013-2015, 2019, 2023 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef SAMTOOLS_H
#define SAMTOOLS_H

#include "htslib/hts_defs.h"
#include "htslib/sam.h"
#include "sam_utils.h"

const char *samtools_version(void);

/* BAM sanitizer options */
#define FIX_POS   2
#define FIX_MQUAL 4
#define FIX_UNMAP 8
#define FIX_CIGAR 16
#define FIX_AUX   32

// default for position sorted data
#define FIX_ON (FIX_MQUAL|FIX_UNMAP|FIX_CIGAR|FIX_AUX)
#define FIX_ALL 255

// Parses a comma-separated list of "pos", "mqual", "unmap", "cigar", and "aux"
// keywords for the bam sanitizer.
int bam_sanitize_options(const char *str);

// Sanitize a BAM record, using FIX_* bit flags as defined above.
// Returns 0 on success,
//        <0 on failure.
int bam_sanitize(sam_hdr_t *h, bam1_t *b, int flags);

#endif
