/*  sam_utils.c -- to hold utility functions and types

    Copyright (C) 2023 Genome Research Ltd.

    Author: Vasudeva Sarma <vasudeva.sarma@sanger.ac.uk>

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

#ifndef SAM_UTIL_H
#define SAM_UTIL_H

#include "htslib/khash.h"
#include "htslib/sam.h"

//this file may contain any utility functions and data types to be shared across

/*below parse_aux_list and aux_exists are moved from sam_view.c to here for common
 *usage at different source files
 */

KHASH_SET_INIT_INT(aux_exists)                      //SET data type to hold aux tags
typedef khash_t(aux_exists) *auxhash_t;

/// parse_aux_list - parses given string for aux tags which are ',' separated
/** @param h - pointer to a SET holding aux tags
 * @param optarg - string having the ',' separated aux tags
 * @param msgheader - string to be used during error output as a header
returns -1 on failure and 0 on success
moved from sam_view.c to here for common usage at different source files
*/
int parse_aux_list(auxhash_t *h, char *optarg, const char *msgheader);


// below utility function declarations moved from samtools.h to here and this header is included in samtools.h

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

#endif  //SAM_UTIL_H


