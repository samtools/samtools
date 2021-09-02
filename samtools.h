/*  samtools.h -- utility routines.

    Copyright (C) 2013-2015, 2019 Genome Research Ltd.

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

const char *samtools_version(void);

#define CHECK_PRINTF(fmt,args) HTS_FORMAT(HTS_PRINTF_FMT, (fmt), (args))

void print_error(const char *subcommand, const char *format, ...) CHECK_PRINTF(2, 3);
void print_error_errno(const char *subcommand, const char *format, ...) CHECK_PRINTF(2, 3);

void check_sam_close(const char *subcmd, samFile *fp, const char *fname, const char *null_fname, int *retp);

/* Utility functions to register an output htsFile/samFile/vcfFile that
 * might be stdout. If FNAME is "-" or NULL, records FP so that print_error()
 * et al can automatically flush it before printing an error message.
 */
void autoflush_if_stdout(htsFile *fp, const char *fname);

/* Call this before closing FP; check_sam_close() does this automatically.
 */
void release_autoflush(htsFile *fp);

/*
 * Utility function to add an index to a file we've opened for write.
 * NB: Call this after writing the header and before writing sequences.
 *
 * The returned index filename should be freed by the caller, but only
 * after sam_idx_save has been called.
 *
 * Returns index filename on success,
 *         NULL on failure.
 */
char *auto_index(htsFile *fp, const char *fn, bam_hdr_t *header);

#endif
