/*  test/merge/test_trans_tbl_init.c -- merge test harness.

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
#include <assert.h>
#include <regex.h>

typedef struct refseq_info {
    const char *name;
    uint32_t    len;
} refseq_info_t;

void dump_header(bam_hdr_t* hdr) {
    printf("->n_targets:(%d)\n", hdr->n_targets);
    int i;
    for (i = 0; i < hdr->n_targets; ++i) {
        printf("->target_name[%d]:(%s)\n",i,hdr->target_name[i]);
        printf("->target_len[%d]:(%d)\n",i,hdr->target_len[i]);
    }

    printf("->text:(");
    fwrite((void*)hdr->text, (size_t) hdr->l_text, 1, stdout);
    printf(")\n");
}

static int populate_merged_header(bam_hdr_t *hdr, merged_header_t *merged_hdr) {
    trans_tbl_t dummy;
    int res;
    res = trans_tbl_init(merged_hdr, hdr, &dummy, 0, 0, 1, NULL);
    trans_tbl_destroy(&dummy);
    return res;
}

/*
 * Populate merged_hdr with data from bam0_header_text and bam0_refseqs.
 * Return bam_hdr_t based on the content in bam1_header_text and bam1_refseqs.
 */

bam_hdr_t * setup_test(const char *bam0_header_text,
                       const refseq_info_t *bam0_refseqs,
                       int32_t bam0_n_refseqs,
                       const char *bam1_header_text,
                       const refseq_info_t *bam1_refseqs,
                       int32_t bam1_n_refseqs,
                       merged_header_t *merged_hdr) {
    bam_hdr_t* bam0 = NULL;
    bam_hdr_t* bam1 = NULL;
    int32_t i;

    bam0 = bam_hdr_init();
    bam0->text = strdup(bam0_header_text);
    if (!bam0->text) goto fail;
    bam0->l_text = strlen(bam0_header_text);
    bam0->n_targets = 1;
    bam0->target_name = (char**)calloc(bam0_n_refseqs, sizeof(char*));
    bam0->target_len = (uint32_t*)calloc(bam0_n_refseqs, sizeof(uint32_t));
    for (i = 0; i < bam0_n_refseqs; i++) {
        bam0->target_name[i] = strdup(bam0_refseqs[i].name);
        if (!bam0->target_name[i]) goto fail;
        bam0->target_len[i] = bam0_refseqs[i].len;
    }

    if (populate_merged_header(bam0, merged_hdr)) goto fail;

    bam1 = bam_hdr_init();
    if (!bam1) goto fail;
    bam1->text = strdup(bam1_header_text);
    if (!bam1->text) goto fail;
    bam1->l_text = strlen(bam1_header_text);
    bam1->n_targets = bam1_n_refseqs;
    bam1->target_name = (char**)calloc(bam1_n_refseqs, sizeof(char*));
    bam1->target_len = (uint32_t*)calloc(bam1_n_refseqs, sizeof(uint32_t));
    for (i = 0; i < bam1_n_refseqs; i++) {
        bam1->target_name[i] = strdup(bam1_refseqs[i].name);
        if (!bam1->target_name[i]) goto fail;
        bam1->target_len[i] = bam1_refseqs[i].len;
    }

    bam_hdr_destroy(bam0);
    return bam1;

 fail:
    bam_hdr_destroy(bam1);
    bam_hdr_destroy(bam0);
    return NULL;
}

#define NELE(x) (sizeof((x)) / sizeof((x)[0]))

static const char init_text[] =
    "@HD\tVN:1.4\tSO:unknown\n"
    "@SQ\tSN:fish\tLN:133\tSP:frog";

static const refseq_info_t init_refs[1] = {
    { "fish", 133 }
};

static const char test_1_trans_text[] =
"@HD\tVN:1.4\tSO:unknown\n"
"@SQ\tSN:fish\tLN:133\n";

static const refseq_info_t test_1_refs[1] = {
    { "fish", 133 }
};

bam_hdr_t * setup_test_1(merged_header_t *merged_hdr) {
    return setup_test(init_text, init_refs, NELE(init_refs),
                      test_1_trans_text, test_1_refs, NELE(test_1_refs),
                      merged_hdr);
}

bool check_test_1(bam_hdr_t* translate, bam_hdr_t* out, trans_tbl_t* tbl) {
    // Check input is unchanged
    if (
        strncmp(test_1_trans_text, translate->text, translate->l_text)
        || translate->l_text != strlen( test_1_trans_text)
        || translate->n_targets != 1
        ) return false;

    // Check output header
    const char out_regex[] =
    "^@HD\tVN:1.4\tSO:unknown\n"
    "@SQ\tSN:fish\tLN:133\tSP:frog\n$";

    regex_t check_regex;
    regcomp(&check_regex, out_regex, REG_EXTENDED|REG_NOSUB);

    if ( regexec(&check_regex, out->text, 0, NULL, 0) != 0 || out->n_targets != 1 ) return false;

    regfree(&check_regex);

    // Check output tbl
    if (tbl[0].n_targets != 1 || tbl[0].tid_trans[0] != 0 || tbl[0].lost_coord_sort) return false;

    return true;
}

static const char test_2_trans_text[] =
"@HD\tVN:1.4\tSO:unknown\n"
"@SQ\tSN:donkey\tLN:133\n"
"@SQ\tSN:fish\tLN:133";

static const refseq_info_t test_2_refs[2] = {
    { "donkey", 133 },
    { "fish",   133 }
};

bam_hdr_t * setup_test_2(merged_header_t *merged_hdr) {
    return setup_test(init_text, init_refs, NELE(init_refs),
                      test_2_trans_text, test_2_refs, NELE(test_2_refs),
                      merged_hdr);
}

bool check_test_2(bam_hdr_t* translate, bam_hdr_t* out, trans_tbl_t* tbl) {
    // Check input is unchanged
    if (
        strncmp(test_2_trans_text, translate->text, translate->l_text)
        || translate->l_text != strlen(test_2_trans_text)
        || translate->n_targets != 2
        ) return false;

    // Check output header
    const char out_regex[] =
    "^@HD\tVN:1.4\tSO:unknown\n"
    "@SQ\tSN:fish\tLN:133\tSP:frog\n"
    "@SQ\tSN:donkey\tLN:133\n$";

    regex_t check_regex;
    regcomp(&check_regex, out_regex, REG_EXTENDED|REG_NOSUB);

    if ( regexec(&check_regex, out->text, 0, NULL, 0) != 0 || out->n_targets != 2 ) return false;

    regfree(&check_regex);

    // Check output tbl
    if (tbl[0].n_targets != 2 || tbl[0].tid_trans[0] != 1 || tbl[0].tid_trans[1] != 0) return false;

    return true;
}

static const char test_3_trans_text[] =
"@HD\tVN:1.4\tSO:unknown\n"
"@SQ\tSN:donkey\tLN:133\n"
"@SQ\tSN:fish\tLN:133\n"
"@RG\tID:fish\tPU:trans\n";

static const refseq_info_t test_3_refs[2] = {
    { "donkey", 133 },
    { "fish",   133 }
};

bam_hdr_t * setup_test_3(merged_header_t *merged_hdr) {
    return setup_test(init_text, init_refs, NELE(init_refs),
                      test_3_trans_text, test_3_refs, NELE(test_3_refs),
                      merged_hdr);
}

bool check_test_3(bam_hdr_t* translate, bam_hdr_t* out, trans_tbl_t* tbl) {
    // Check input is unchanged
    if (
        strncmp(test_3_trans_text, translate->text, translate->l_text)
        || translate->l_text != strlen(test_3_trans_text)
        || translate->n_targets != 2
        ) return false;
    return true;
}

static const char test_4_trans_text[] =
"@HD\tVN:1.4\tSO:unknown\n"
"@SQ\tSN:donkey\tLN:133\n"
"@SQ\tSN:fish\tLN:133\n"
"@RG\tID:fish\tPU:trans\n";

static const refseq_info_t test_4_refs[2] = {
    { "donkey", 133 },
    { "fish",   133 }
};

bam_hdr_t * setup_test_4(merged_header_t *merged_hdr) {
    const char* t4_init_text =
        "@HD\tVN:1.4\tSO:unknown\n"
        "@SQ\tSN:fish\tLN:133\tSP:frog\n"
        "@RG\tID:fish\tPU:out\n";

    return setup_test(t4_init_text, init_refs, NELE(init_refs),
                      test_4_trans_text, test_4_refs, NELE(test_4_refs),
                      merged_hdr);
}

bool check_test_4(bam_hdr_t* translate, bam_hdr_t* out, trans_tbl_t* tbl) {
    // Check input is unchanged
    if (
        strncmp(test_4_trans_text, translate->text, translate->l_text)
        || translate->l_text != strlen(test_4_trans_text)
        || translate->n_targets != 2
        ) return false;
    return true;
}

static const char test_5_trans_text[] =
"@HD\tVN:1.4\tSO:unknown\n"
"@SQ\tSN:donkey\tLN:133\n"
"@SQ\tSN:fish\tLN:133\n"
"@RG\tID:fish\tPU:trans\n"
"@PG\tXX:dummy\tID:fish\tDS:trans\n"
"@PG\tPP:fish\tID:hook\tDS:trans\n";

static const refseq_info_t test_5_refs[2] = {
    { "donkey", 133 },
    { "fish",   133 }
};

bam_hdr_t * setup_test_5(merged_header_t *merged_hdr) {
    const char* t5_init_text =
        "@HD\tVN:1.4\tSO:unknown\n"
        "@SQ\tSN:fish\tLN:133\tSP:frog\n"
        "@RG\tID:fish\tPU:out\n"
        "@PG\tXX:dummyx\tID:fish\tDS:out\n"
        "@PG\tPP:fish\tID:hook\tDS:out\n";

    return setup_test(t5_init_text, init_refs, NELE(init_refs),
                      test_5_trans_text, test_5_refs, NELE(test_5_refs),
                      merged_hdr);
}

bool check_test_5(bam_hdr_t* translate, bam_hdr_t* out, trans_tbl_t* tbl) {
    // Check input is unchanged
    if (
        strncmp(test_5_trans_text, translate->text, translate->l_text)
        || translate->l_text != strlen(test_5_trans_text)
        || translate->n_targets != 2
        ) return false;
    return true;
}

static const char test_6_trans_text[] =
"@HD\tVN:1.4\tSO:unknown\n"
"@SQ\tSN:donkey\tLN:133\n"
"@SQ\tSN:fish\tLN:133\n"
"@RG\tID:fish\tPU:trans\n"
"@PG\tXX:dummy\tID:fish\tDS:trans\n"
"@PG\tPP:fish\tID:hook\tDS:trans\n";

static const refseq_info_t test_6_refs[2] = {
    { "donkey", 133 },
    { "fish",   133 }
};

bam_hdr_t * setup_test_6(merged_header_t *merged_hdr) {
    return setup_test(init_text, init_refs, NELE(init_refs),
                      test_6_trans_text, test_6_refs, NELE(test_6_refs),
                      merged_hdr);
}

bool check_test_6(bam_hdr_t* translate, bam_hdr_t* out, trans_tbl_t* tbl) {
    // Check input is unchanged
    if (
        strncmp(test_6_trans_text, translate->text, translate->l_text)
        || translate->l_text != strlen(test_5_trans_text)
        || translate->n_targets != 2
        ) return false;
    return true;
}

int main(int argc, char**argv)
{
    const int NUM_TESTS = 6;
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

    // Set the seed to a fixed value so that calls to lrand48 within functions return predictable values
    const long GIMMICK_SEED = 0x1234330e;
    srand48(GIMMICK_SEED);

    bam_hdr_t* out;
    bam_hdr_t* translate;

    if (verbose) printf("BEGIN test 1\n");
    // setup
    trans_tbl_t tbl_1;
    merged_header_t *merged_hdr = init_merged_header();
    translate = setup_test_1(merged_hdr);
    assert(translate);
    // test
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
    }
    if (verbose) printf("RUN test 1\n");
    trans_tbl_init(merged_hdr, translate, &tbl_1, false, false, true, NULL);
    out = finish_merged_header(merged_hdr);
    free_merged_header(merged_hdr);
    if (verbose) printf("END RUN test 1\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_1(translate, out, &tbl_1)) {
        if (verbose) printf("Test 1 : PASS\n");
        ++success;
    } else {
        if (verbose) printf("Test 1 : FAIL\n");
        fprintf(stderr, "Test 1 : FAIL\n");
        ++failure;
    }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_1);
    if (verbose) printf("END test 1\n");

    // test
    if (verbose) printf("BEGIN test 2\n");
    // reinit
    trans_tbl_t tbl_2;

    merged_hdr = init_merged_header();
    translate = setup_test_2(merged_hdr);
    assert(translate);
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
    }
    if (verbose) printf("RUN test 2\n");
    trans_tbl_init(merged_hdr, translate, &tbl_2, false, false, true, NULL);
    out = finish_merged_header(merged_hdr);
    free_merged_header(merged_hdr);
    if (verbose) printf("END RUN test 2\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_2(translate, out, &tbl_2)) {
        if (verbose) printf("Test 2 : PASS\n");
        ++success;
    } else {
        if (verbose) printf("Test 2 : FAIL\n");
        fprintf(stderr, "Test 2 : FAIL\n");
        ++failure;
    }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_2);
    if (verbose) printf("END test 2\n");

    // test
    if (verbose) printf("BEGIN test 3\n");
    // reinit
    trans_tbl_t tbl_3;
    merged_hdr = init_merged_header();
    translate = setup_test_3(merged_hdr);
    assert(translate);
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
     }
    if (verbose) printf("RUN test 3\n");
    trans_tbl_init(merged_hdr, translate, &tbl_3, false, false, true, NULL);
    out = finish_merged_header(merged_hdr);
    free_merged_header(merged_hdr);
    if (verbose) printf("END RUN test 3\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_3(translate, out, &tbl_3)) {
        if (verbose) printf("Test 3 : PASS\n");
        ++success;
    } else {
        if (verbose) printf("Test 3 : FAIL\n");
        fprintf(stderr, "Test 3 : FAIL\n");
        ++failure;
    }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_3);
    if (verbose) printf("END test 3\n");

    // test
    if (verbose) printf("BEGIN test 4\n");
    // reinit
    trans_tbl_t tbl_4;
    merged_hdr = init_merged_header();
    translate = setup_test_4(merged_hdr);
    assert(translate);
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
    }
    if (verbose) printf("RUN test 4\n");
    trans_tbl_init(merged_hdr, translate, &tbl_4, false, false, true, NULL);
    out = finish_merged_header(merged_hdr);
    free_merged_header(merged_hdr);
    if (verbose) printf("END RUN test 4\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_4(translate, out, &tbl_4)) {
        if (verbose) printf("Test 4 : PASS\n");
        ++success;
    } else {
        if (verbose) printf("Test 4 : FAIL\n");
        fprintf(stderr, "Test 4 : FAIL\n");
        ++failure;
    }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_4);
    if (verbose) printf("END test 4\n");

    // test
    if (verbose) printf("BEGIN test 5\n");
    // reinit
    trans_tbl_t tbl_5;
    merged_hdr = init_merged_header();
    translate = setup_test_5(merged_hdr);
    assert(translate);
    if (verbose > 1) {

        printf("translate\n");
        dump_header(translate);
    }
    if (verbose) printf("RUN test 5\n");
    trans_tbl_init(merged_hdr, translate, &tbl_5, false, false, true, NULL);
    out = finish_merged_header(merged_hdr);
    free_merged_header(merged_hdr);
    if (verbose) printf("END RUN test 5\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_5(translate, out, &tbl_5)) {
        if (verbose) printf("Test 5 : PASS\n");
        ++success;
    } else {
        if (verbose) printf("Test 5 : FAIL\n");
        fprintf(stderr, "Test 5 : FAIL\n");
        ++failure;
    }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_5);
    if (verbose) printf("END test 5\n");

    // test
    if (verbose) printf("BEGIN test 6\n");
    // reinit
    trans_tbl_t tbl_6;
    merged_hdr = init_merged_header();
    translate = setup_test_6(merged_hdr);
    assert(translate);
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
    }
    if (verbose) printf("RUN test 6\n");
    trans_tbl_init(merged_hdr, translate, &tbl_6, false, false, true, "filename");
    out = finish_merged_header(merged_hdr);
    free_merged_header(merged_hdr);
    if (verbose) printf("END RUN test 6\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_6(translate, out, &tbl_6)) {
        if (verbose) printf("Test 6 : PASS\n");
        ++success;
    } else {
        if (verbose) printf("Test 6 : FAIL\n");
        fprintf(stderr, "Test 6 : FAIL\n");
        ++failure;
    }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_6);
    if (verbose) printf("END test 6\n");

    if (success == NUM_TESTS) {
        return 0;
    } else {
        fprintf(stderr, "%d failures %d successes\n", failure, success);
        return 1;
    }
}
