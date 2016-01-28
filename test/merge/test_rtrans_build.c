/*  test/merge/test_rtrans_build.c -- header translation test harness.

    Copyright (C) 2013, 2014 Genome Research Ltd.

    Author: Martin O. Pollard <mp15@sanger.ac.uk>

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

#include "../../bam_sort.c"

void dump_rtrans(int* rtrans, int n, int n_targets) {
    printf("->n_targets:(%d)\n", n_targets);
    int i, j;
    for (i = 0; i < n; ++i) {
        fprintf(stderr, "%d",rtrans[i*n_targets+0]);
        for (j = 1; j < n_targets; ++j)
            fprintf(stderr, "\t%d",rtrans[i*n_targets+j]);
        fprintf(stderr, "\n");
    }
}

void setup_test_1(trans_tbl_t* tbl) {
    tbl[0].n_targets = 2;
    tbl[0].tid_trans = calloc(sizeof(int), 2);
    tbl[0].tid_trans[0] = 0;
    tbl[0].tid_trans[1] = 1;
    tbl[0].rg_trans = kh_init(c2c);
    tbl[0].pg_trans = kh_init(c2c);

    tbl[1].n_targets = 2;
    tbl[1].tid_trans = calloc(sizeof(int), 2);
    tbl[1].tid_trans[0] = 1;
    tbl[1].tid_trans[1] = 2;
    tbl[1].rg_trans = kh_init(c2c);
    tbl[1].pg_trans = kh_init(c2c);
}

bool check_test_1(trans_tbl_t* tbl, int* rtrans) {
    // Check input is unchanged

    // Check output

    return true;
}


int main(int argc, char**argv)
{
    const int NUM_TESTS = 1;
    int verbose = 0;
    int success = 0;
    int failure = 0;
    int getopt_char;
    while ((getopt_char = getopt(argc, argv, "v")) != -1) {
        switch (getopt_char) {
            case 'v':
                ++verbose;
                break;
            default:
                break;
        }
    }
    const long GIMMICK_SEED = 0x1234330e;
    srand48(GIMMICK_SEED);

    if (verbose) printf("BEGIN test 1\n");
    // setup
    trans_tbl_t tbl_1[2];
    int n_targets_1 = 3;
    int n_1 = 2;
    int* rtrans_1 = NULL;
    setup_test_1(&tbl_1[0]);
    // test
    if (verbose > 1) {
        // dump_trans_tid
    }
    if (verbose) printf("RUN test 1\n");
    rtrans_1 = rtrans_build(n_1, n_targets_1, &tbl_1[0]);
    if (verbose) printf("END RUN test 1\n");
    if (verbose > 1) {
        printf("rtrans\n");
        dump_rtrans(rtrans_1, n_1, n_targets_1);
    }
    if (check_test_1(&tbl_1[0], rtrans_1)) {
        ++success;
    } else {
        ++failure;
        if (verbose) printf("FAIL test 1\n");
    }
    // teardown
    trans_tbl_destroy(&tbl_1[0]);
    trans_tbl_destroy(&tbl_1[1]);
    free(rtrans_1);
    if (verbose) printf("END test 1\n");

    if (success == NUM_TESTS) {
        return 0;
    } else {
        fprintf(stderr, "%d failures %d successes\n", failure, success);
        return 1;
    }
}
