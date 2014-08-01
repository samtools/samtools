/*  stats_isize.h -- generalised insert size calculation for samtools stats.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Nicholas Clarke <nc6@sanger.ac.uk>

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

#include <htslib/khash.h>
#include <stdint.h>

typedef struct
{
    int total;
    uint64_t *isize_inward, *isize_outward, *isize_other;
}
isize_dense_data_t;

typedef struct
{
    uint64_t isize_inward, isize_outward, isize_other;
}
isize_sparse_record_t;

KHASH_MAP_INIT_INT(m32, isize_sparse_record_t *)

typedef struct
{
    int max;
    khash_t(m32) *array;
}
isize_sparse_data_t;

typedef union {
    isize_sparse_data_t *sparse;
    isize_dense_data_t *dense;
} isize_data_t;

// Insert size structure
typedef struct
{
    isize_data_t data;

    // Maximum
    int (*nitems)(isize_data_t);

    // Fetch the number of inserts of a given size
    uint64_t (*inward)(isize_data_t, int);
    uint64_t (*outward)(isize_data_t, int);
    uint64_t (*other)(isize_data_t, int);

    // Set the number of inserts of a given size
    void (*set_inward)(isize_data_t, int, uint64_t);
    void (*set_outward)(isize_data_t, int, uint64_t);
    void (*set_other)(isize_data_t, int, uint64_t);

    // Increment the number of inserts of a given size
    void (*inc_inward)(isize_data_t, int);
    void (*inc_outward)(isize_data_t, int);
    void (*inc_other)(isize_data_t, int);

    // Free this structure
    void (*isize_free)(isize_data_t);
}
isize_t;

isize_t *init_isize_t(int bound);
