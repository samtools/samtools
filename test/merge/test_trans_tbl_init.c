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

#include "../../bam_sort.c"

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

static const char test_1_trans_text[] =
"@HD\tVN:1.4\tSO:unknown\n"
"@SQ\tSN:fish\tLN:133\n";

void setup_test_1(bam_hdr_t** translate_in, bam_hdr_t** out_in) {
    bam_hdr_t* out;
    bam_hdr_t* translate;

    translate = bam_hdr_init();
    translate->text = strdup(test_1_trans_text);
    translate->l_text = strlen(test_1_trans_text);
    translate->n_targets = 1;
    translate->target_name = (char**)calloc(1, sizeof(char*));
    translate->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
    translate->target_name[0] = strdup("fish");
    translate->target_len[0] = 133;
    out = bam_hdr_init();
    const char out_text[] =
        "@HD\tVN:1.4\tSO:unknown\n"
        "@SQ\tSN:fish\tLN:133\tSP:frog";
    out->text = strdup(out_text);
    out->l_text = strlen(out_text);
    out->n_targets = 1;
    out->target_name = (char**)calloc(1, sizeof(char*));
    out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
    out->target_name[0] = strdup("fish");
    out->target_len[0] = 133;

    *translate_in = translate;
    *out_in = out;
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

void setup_test_2(bam_hdr_t** translate_in, bam_hdr_t** out_in) {
    bam_hdr_t* out;
    bam_hdr_t* translate;

    translate = bam_hdr_init();
    translate->text = strdup(test_2_trans_text);
    translate->l_text = strlen(test_2_trans_text);
    translate->n_targets = 2;
    translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
    translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
    translate->target_name[0] = strdup("donkey");
    translate->target_len[0] = 133;
    translate->target_name[1] = strdup("fish");
    translate->target_len[1] = 133;
    out = bam_hdr_init();
    const char* out_text =
        "@HD\tVN:1.4\tSO:unknown\n"
        "@SQ\tSN:fish\tLN:133\tSP:frog";
    out->text = strdup(out_text);
    out->l_text = strlen(out_text);
    out->n_targets = 1;
    out->target_name = (char**)calloc(1, sizeof(char*));
    out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
    out->target_name[0] = strdup("fish");
    out->target_len[0] = 133;

    *translate_in = translate;
    *out_in = out;
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

void setup_test_3(bam_hdr_t** translate_in, bam_hdr_t** out_in) {
    bam_hdr_t* out;
    bam_hdr_t* translate;

    translate = bam_hdr_init();
    translate->text = strdup(test_3_trans_text);
    translate->l_text = strlen(test_3_trans_text);
    translate->n_targets = 2;
    translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
    translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
    translate->target_name[0] = strdup("donkey");
    translate->target_len[0] = 133;
    translate->target_name[1] = strdup("fish");
    translate->target_len[1] = 133;
    out = bam_hdr_init();
    const char* out_text =
        "@HD\tVN:1.4\tSO:unknown\n"
        "@SQ\tSN:fish\tLN:133\tSP:frog";
    out->text = strdup(out_text);
    out->l_text = strlen(out_text);
    out->n_targets = 1;
    out->target_name = (char**)calloc(1, sizeof(char*));
    out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
    out->target_name[0] = strdup("fish");
    out->target_len[0] = 133;

    *translate_in = translate;
    *out_in = out;
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

void setup_test_4(bam_hdr_t** translate_in, bam_hdr_t** out_in) {
    bam_hdr_t* out;
    bam_hdr_t* translate;

    translate = bam_hdr_init();
    translate->text = strdup(test_4_trans_text);
    translate->l_text = strlen(test_4_trans_text);
    translate->n_targets = 2;
    translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
    translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
    translate->target_name[0] = strdup("donkey");
    translate->target_len[0] = 133;
    translate->target_name[1] = strdup("fish");
    translate->target_len[1] = 133;
    out = bam_hdr_init();
    const char* out_text =
        "@HD\tVN:1.4\tSO:unknown\n"
        "@SQ\tSN:fish\tLN:133\tSP:frog\n"
        "@RG\tID:fish\tPU:out\n";
    out->text = strdup(out_text);
    out->l_text = strlen(out_text);
    out->n_targets = 1;
    out->target_name = (char**)calloc(1, sizeof(char*));
    out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
    out->target_name[0] = strdup("fish");
    out->target_len[0] = 133;

    *translate_in = translate;
    *out_in = out;
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

void setup_test_5(bam_hdr_t** translate_in, bam_hdr_t** out_in) {
    bam_hdr_t* out;
    bam_hdr_t* translate;

    translate = bam_hdr_init();
    translate->text = strdup(test_5_trans_text);
    translate->l_text = strlen(test_5_trans_text);
    translate->n_targets = 2;
    translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
    translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
    translate->target_name[0] = strdup("donkey");
    translate->target_len[0] = 133;
    translate->target_name[1] = strdup("fish");
    translate->target_len[1] = 133;
    out = bam_hdr_init();
    const char* out_text =
        "@HD\tVN:1.4\tSO:unknown\n"
        "@SQ\tSN:fish\tLN:133\tSP:frog\n"
        "@RG\tID:fish\tPU:out\n"
        "@PG\tXX:dummyx\tID:fish\tDS:out\n"
        "@PG\tPP:fish\tID:hook\tDS:out\n";
    out->text = strdup(out_text);
    out->l_text = strlen(out_text);
    out->n_targets = 1;
    out->target_name = (char**)calloc(1, sizeof(char*));
    out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
    out->target_name[0] = strdup("fish");
    out->target_len[0] = 133;

    *translate_in = translate;
    *out_in = out;
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


int main(int argc, char**argv)
{
    const int NUM_TESTS = 5;
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
    setup_test_1(&translate,&out);
    // test
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (verbose) printf("RUN test 1\n");
    trans_tbl_init(out, translate, &tbl_1, false, false);
    if (verbose) printf("END RUN test 1\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_1(translate, out, &tbl_1)) { ++success; } else { ++failure; }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_1);
    if (verbose) printf("END test 1\n");

    // test
    if (verbose) printf("BEGIN test 2\n");
    // reinit
    trans_tbl_t tbl_2;
    setup_test_2(&translate,&out);
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (verbose) printf("RUN test 2\n");
    trans_tbl_init(out, translate, &tbl_2, false, false);
    if (verbose) printf("END RUN test 2\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_2(translate, out, &tbl_2)) { ++success; } else { ++failure; }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_2);
    if (verbose) printf("END test 2\n");

    // test
    if (verbose) printf("BEGIN test 3\n");
    // reinit
    trans_tbl_t tbl_3;
    setup_test_3(&translate,&out);
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (verbose) printf("RUN test 3\n");
    trans_tbl_init(out, translate, &tbl_3, false, false);
    if (verbose) printf("END RUN test 3\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_3(translate, out, &tbl_3)) { ++success; } else { ++failure; }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_3);
    if (verbose) printf("END test 3\n");

    // test
    if (verbose) printf("BEGIN test 4\n");
    // reinit
    trans_tbl_t tbl_4;
    setup_test_4(&translate,&out);
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (verbose) printf("RUN test 4\n");
    trans_tbl_init(out, translate, &tbl_4, false, false);
    if (verbose) printf("END RUN test 4\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_4(translate, out, &tbl_4)) { ++success; } else { ++failure; }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_4);
    if (verbose) printf("END test 4\n");

    // test
    if (verbose) printf("BEGIN test 5\n");
    // reinit
    trans_tbl_t tbl_5;
    setup_test_5(&translate,&out);
    if (verbose > 1) {

        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (verbose) printf("RUN test 5\n");
    trans_tbl_init(out, translate, &tbl_5, false, false);
    if (verbose) printf("END RUN test 5\n");
    if (verbose > 1) {
        printf("translate\n");
        dump_header(translate);
        printf("out\n");
        dump_header(out);
    }
    if (check_test_5(translate, out, &tbl_5)) { ++success; } else { ++failure; }
    // teardown
    bam_hdr_destroy(translate);
    bam_hdr_destroy(out);
    trans_tbl_destroy(&tbl_5);
    if (verbose) printf("END test 5\n");

    if (success == NUM_TESTS) {
        return 0;
    } else {
        fprintf(stderr, "%d failures %d successes\n", failure, success);
        return 1;
    }
}
