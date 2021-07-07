/*  bam_ampliconclip.h -- shared functions between amplicon clip/stats

    Copyright (C) 2020-2021 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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

#ifndef BAM_AMPLICONCLIP_H
#define BAM_AMPLICONCLIP_H

#include "htslib/khash.h"

typedef struct {
    int64_t left;
    int64_t right;
    int rev;
} bed_entry_t;

typedef struct {
    bed_entry_t *bp;
    int64_t longest;
    int length;
    int size;
} bed_entry_list_t;

KHASH_MAP_INIT_STR(bed_list_hash, bed_entry_list_t);

#define BED_LIST_INIT {NULL, 0, 0, 0, {0}}


int load_bed_file_multi_ref(char *infile, int get_strand,
                        int sort_by_pos, khash_t(bed_list_hash) *bed_lists);

void destroy_bed_hash(khash_t(bed_list_hash) *hash);


#endif /* BAM_AMPLICONCLIP_H */
