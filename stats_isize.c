/*  stats_isize.c -- generalised insert size calculation for samtools stats.

    Copyright (C) 2014, 2018 Genome Research Ltd.

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

#include <config.h>

#include <stdio.h>
#include "stats_isize.h"
#include <htslib/khash.h>

typedef enum {IN,OUT,OTHER} isize_insert_t;

static int max(int a, int b) {
    if (a < b) {
        return b;
    } else {
        return a;
    }
}

static isize_sparse_record_t * sparse_get_f(isize_data_t data, int at) {
    isize_sparse_data_t *a = data.sparse;
    khash_t(m32) *h = a->array;

    khint_t k = kh_get(m32, h, at);
    if (k != kh_end(h)) {
        return kh_value(h, k);
    } else {
        return NULL;
    }
}

static uint64_t sparse_in_f(isize_data_t data, int at) {
    isize_sparse_record_t* a = sparse_get_f(data, at);
    if (a != NULL) {
        return a->isize_inward;
    } else {
        return 0;
    }
}
static uint64_t sparse_out_f(isize_data_t data, int at) {
    isize_sparse_record_t* a = sparse_get_f(data, at);
    if (a != NULL) {
        return a->isize_outward;
    } else {
        return 0;
    }
}
static uint64_t sparse_other_f(isize_data_t data, int at) {
    isize_sparse_record_t* a = sparse_get_f(data, at);
    if (a != NULL) {
        return a->isize_other;
    } else {
        return 0;
    }
}

static void sparse_set_f(isize_data_t data, int at, isize_insert_t field, uint64_t value) {
    isize_sparse_data_t *a = data.sparse;
    khash_t(m32) *h = a->array;

    khint_t k = kh_get(m32, h, at);
    isize_sparse_record_t *rec;
    if (k != kh_end(h)) {
        rec = kh_value(h, k);
    } else if (value != 0) {
        rec = malloc(sizeof(isize_sparse_record_t));
        if (rec != NULL) {
            rec->isize_inward = 0;
            rec->isize_outward = 0;
            rec->isize_other = 0;
            int stupid = 0;
            khint_t it = kh_put(m32, h, at, & stupid);
            kh_value(h, it) = rec;
            a->max = max(at, a->max);
        } else {
            fprintf(stderr, "%s\n", "Failed to allocate memory for isize_sparse_record_t");
            exit(11);
        }
    } else {
        return;
    }
    if (field == IN) {
        rec->isize_inward = value;
    } else if (field == OUT) {
        rec->isize_outward = value;
    } else {
        rec->isize_other = value;
    }

}

static void sparse_set_in_f(isize_data_t data, int at, uint64_t value) { sparse_set_f(data, at, IN, value); }
static void sparse_set_out_f(isize_data_t data, int at, uint64_t value) { sparse_set_f(data, at, OUT, value); }
static void sparse_set_other_f(isize_data_t data, int at, uint64_t value) { sparse_set_f(data, at, OTHER, value); }

static void sparse_inc_in_f(isize_data_t data, int at) { sparse_set_in_f(data, at, sparse_in_f(data, at) + 1); }
static void sparse_inc_out_f(isize_data_t data, int at) { sparse_set_out_f(data, at, sparse_out_f(data, at) + 1); }
static void sparse_inc_other_f(isize_data_t data, int at) { sparse_set_other_f(data, at, sparse_other_f(data, at) + 1); }

static void sparse_isize_free(isize_data_t data) {
    isize_sparse_data_t *a = data.sparse;
    khint_t k;
    for (k = 0; k < kh_end(a->array); ++k)
        if (kh_exist(a->array, k)) free(kh_val(a->array, k));
    kh_destroy(m32, a->array);
    free(a);
}

static int sparse_nitems(isize_data_t data) {
    isize_sparse_data_t *a = data.sparse;
    return a->max + 1;
}

static uint64_t dense_in_f(isize_data_t data, int at) { return data.dense->isize_inward[at]; }
static uint64_t dense_out_f(isize_data_t data, int at) { return data.dense->isize_outward[at]; }
static uint64_t dense_other_f(isize_data_t data, int at) { return data.dense->isize_other[at]; }

static void dense_set_in_f(isize_data_t data, int at, uint64_t value) { data.dense->isize_inward[at] = value; }
static void dense_set_out_f(isize_data_t data, int at, uint64_t value) { data.dense->isize_outward[at] = value; }
static void dense_set_other_f(isize_data_t data, int at, uint64_t value) { data.dense->isize_other[at] = value; }

static void dense_inc_in_f(isize_data_t data, int at) { data.dense->isize_inward[at] += 1; }
static void dense_inc_out_f(isize_data_t data, int at) { data.dense->isize_outward[at] += 1; }
static void dense_inc_other_f(isize_data_t data, int at) { data.dense->isize_other[at] += 1; }

static void dense_isize_free(isize_data_t data) {
    isize_dense_data_t *a = data.dense;
    free(a->isize_inward);
    free(a->isize_outward);
    free(a->isize_other);
    free(a);
}

static int dense_nitems(isize_data_t data) {
    isize_dense_data_t *a = data.dense;
    return a->total;
}

// Construct a relevant isize_t given the bound.
isize_t *init_isize_t(int bound) {
    if (bound <= 0) {
        // Use sparse data structure.
        isize_sparse_data_t *data = (isize_sparse_data_t *) malloc(sizeof(isize_sparse_data_t));
        if (!data)
            return NULL;

        // Initialise
        data->max = 0;
        data->array = kh_init(m32);
        if (!data->array) {
            free(data);
            return NULL;
        }

        isize_t *isize = (isize_t *)malloc(sizeof(isize_t));
        if (!isize) {
            kh_destroy(m32, data->array);
            free(data);
            return NULL;
        }

        isize->data.sparse = data;
        isize->nitems = & sparse_nitems;

        isize->inward = & sparse_in_f;
        isize->outward = & sparse_out_f;
        isize->other = & sparse_other_f;

        isize->set_inward = & sparse_set_in_f;
        isize->set_outward = & sparse_set_out_f;
        isize->set_other = & sparse_set_other_f;

        isize->inc_inward = & sparse_inc_in_f;
        isize->inc_outward = & sparse_inc_out_f;
        isize->inc_other = & sparse_inc_other_f;

        isize->isize_free = & sparse_isize_free;

        return isize;
    } else {
        uint64_t* in = calloc(bound,sizeof(uint64_t));
        uint64_t* out = calloc(bound,sizeof(uint64_t));
        uint64_t* other = calloc(bound,sizeof(uint64_t));
        isize_dense_data_t *rec = (isize_dense_data_t *)malloc(sizeof(isize_dense_data_t));
        isize_t *isize = (isize_t *)malloc(sizeof(isize_t));
        if (!in || !out || !other || !rec || !isize) {
            free(in);
            free(out);
            free(other);
            free(rec);
            free(isize);
            return NULL;
        }
        rec->isize_inward = in;
        rec->isize_outward = out;
        rec->isize_other = other;
        rec->total=bound;

        isize->data.dense = rec;
        isize->nitems = & dense_nitems;

        isize->inward = & dense_in_f;
        isize->outward = & dense_out_f;
        isize->other = & dense_other_f;

        isize->set_inward = & dense_set_in_f;
        isize->set_outward = & dense_set_out_f;
        isize->set_other = & dense_set_other_f;

        isize->inc_inward = & dense_inc_in_f;
        isize->inc_outward = & dense_inc_out_f;
        isize->inc_other = & dense_inc_other_f;

        isize->isize_free = & dense_isize_free;

        return isize;
    }
}
