/*  bedidx.h -- BED file indexing header file.

    Copyright (C) 2017 Genome Research Ltd.

    Author: Valeriu Ohan <vo2@sanger.ac.uk>

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

#ifndef BEDIDX_H
#define BEDIDX_H

#include "htslib/hts.h"

#define LIDX_SHIFT 13
#define ALL 0
#define FILTERED 1

#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )

/*
 * Bed file reading.  The main functions are bed_read and bed_destroy.
 * The void* return value from bed_read is the same "void *reg_hash" used in
 * bed_overlap, bed_hash_regions, bed_get and bed_reglist.
 *
 * Due to the data structures used, the order of chromosomes in the BED file
 * is not preserved and iterators will return chromosomes in hash key order.
 */
void *bed_read(const char *fn);
void bed_destroy(void *_h);

/*
 * bed_reglist will turn the bed file into a list of regions structured as
 * an hts_reglist_t: an array with one entry per chromosome, with each
 * element containing a further array of regions within that chromosome.
 *
 * Note however this function does not take a SAM or VCF header, so the
 * hts_reglist_t "tid" member will be invalid.  Users are expected to perform
 * this conversion themselves from the "chr" string member instead.  (See
 * hts_itr_regions for example code which does this.)
 *
 * Use filtered = ALL (0) unless you have previously applied bed_hash_regions
 * with reg_hash != NULL and *op == 0, in which case FILTERED may be used to
 * perform an intersection of bed file and an array of region strings.
 *
 * Free using hts_reglist_free()
 */
hts_reglist_t *bed_reglist(void *reg_hash, int filter, int *count_regs);

int bed_overlap(const void *_h, const char *chr, hts_pos_t beg, hts_pos_t end);
void *bed_hash_regions(void *reg_hash, char **regs, int first, int last, int *op);
void bed_unify(void *_h);

#endif
