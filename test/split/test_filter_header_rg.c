/*  test/split/test_filter_header_rg.c -- split test cases.

    Copyright (C) 2014-2016, 2018, 2019 Genome Research Ltd.

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

#include "../test.h"
#include <unistd.h>
#include <stdbool.h>
#include "samtools.h"
#include <string.h>
#include <stdlib.h>
#include "htslib/kstring.h"

int line_cmp(const void *av, const void *bv) {
    const char *a = *(const char **) av;
    const char *b = *(const char **) bv;
    size_t al = strcspn(a, "\n");
    size_t bl = strcspn(b, "\n");
    size_t min = al < bl ? al : bl;
    int m = memcmp(a, b, min);
    if (m != 0) return m;
    if (al < bl) return -1;
    return al == bl ? 0 : 1;
}

bool hdrcmp(const char *hdr1, const char *hdr2) {
    size_t nl1, nl2, count1 = 0, count2 = 0, i;
    const char *l;
    const char **lines1, **lines2;
    int res = 0;

    // First line should be @HD
    if (strncmp(hdr1, "@HD\t", 4) != 0) return false;
    if (strncmp(hdr2, "@HD\t", 4) != 0) return false;
    nl1 = strcspn(hdr1, "\n");
    nl2 = strcspn(hdr2, "\n");
    if (nl1 != nl2 || memcmp(hdr1, hdr2, nl1) != 0) return false;

    // Count lines.
    for (l = hdr1 + nl1; *l != '\0'; l += strcspn(l, "\n")) ++l, ++count1;
    for (l = hdr2 + nl2; *l != '\0'; l += strcspn(l, "\n")) ++l, ++count2;
    if (count1 != count2) return false;

    lines1 = malloc(count1 * sizeof(*lines1));
    if (!lines1) return false;
    lines2 = malloc(count2 * sizeof(*lines2));
    if (!lines2) { free(lines1); return false; }

    for (i = 0, l = hdr1 + nl1; *l != '\0'; l += strcspn(l, "\n"))
        lines1[i++] = ++l;
    for (i = 0, l = hdr2 + nl2; *l != '\0'; l += strcspn(l, "\n"))
        lines2[i++] = ++l;

    qsort(lines1, count1, sizeof(*lines1), line_cmp);
    qsort(lines2, count2, sizeof(*lines2), line_cmp);

    for (i = 0; i < count1; i++) {
        res = line_cmp(&lines1[i], &lines2[i]);
        if (res != 0) break;
    }

    free(lines1);
    free(lines2);

    return res?false:true;
}

void setup_test_1(sam_hdr_t** hdr_in)
{
    *hdr_in = sam_hdr_init();
    const char *test1 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\tLN:1\n"
    "@RG\tID:fish\n";
    sam_hdr_add_lines(*hdr_in, test1, 0);
}

bool check_test_1(sam_hdr_t* hdr) {
    const char *test1_res =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\tLN:1\n"
    "@PG\tID:samtools\tPN:samtools\tVN:x.y.test\tCL:test_filter_header_rg foo bar baz\n";

    return hdrcmp(sam_hdr_str(hdr), test1_res);
}

void setup_test_2(sam_hdr_t** hdr_in)
{
    *hdr_in = sam_hdr_init();
    const char *test2 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\tLN:1\n"
    "@RG\tID:fish\n";
    sam_hdr_add_lines(*hdr_in, test2, 0);
}

bool check_test_2(sam_hdr_t* hdr) {
    const char *test2_res =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\tLN:1\n"
    "@RG\tID:fish\n"
    "@PG\tID:samtools\tPN:samtools\tVN:x.y.test\tCL:test_filter_header_rg foo bar baz\n";

    return hdrcmp(sam_hdr_str(hdr), test2_res);
}

void setup_test_3(sam_hdr_t** hdr_in)
{
    *hdr_in = sam_hdr_init();
    const char *test3 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\tLN:1\n"
    "@RG\tID:fish1\n"
    "@RG\tID:fish2\n"
    "@RG\tID:fish3\n"
    "@RG\tID:fish4\n";
    sam_hdr_add_lines(*hdr_in, test3, 0);
}

bool check_test_3(sam_hdr_t* hdr) {
    const char *test3_res =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\tLN:1\n"
    "@PG\tID:samtools\tPN:samtools\tVN:x.y.test\tCL:test_filter_header_rg foo bar baz\n";

    return hdrcmp(sam_hdr_str(hdr), test3_res);
}

int main(int argc, char *argv[])
{
    // test state
    const int NUM_TESTS = 3;
    int verbose = 0;
    int success = 0;
    int failure = 0;

    int getopt_char;
    char *test_argv[] = { "test_filter_header_rg", "foo\tbar", "baz" };
    char *arg_list = stringify_argv(3, test_argv);
    while ((getopt_char = getopt(argc, argv, "v")) != -1) {
        switch (getopt_char) {
            case 'v':
                ++verbose;
                break;
            default:
                printf(
                       "usage: test_filter_header_rg [-v]\n\n"
                       " -v verbose output\n"
                       );
                break;
        }
    }


    // Setup stderr redirect
    kstring_t res = { 0, 0, NULL };
    int orig_stderr = dup(STDERR_FILENO); // Save stderr
    int redirected_stderr;
    char* tempfname = (optind < argc)? argv[optind] : "test_count_rg.tmp";
    FILE* check = NULL;

    // setup
    if (verbose) printf("BEGIN test 1\n");  // test eliminating a tag that isn't there
    sam_hdr_t* hdr1;
    const char* id_to_keep_1 = "1#2.3";
    setup_test_1(&hdr1);
    if (verbose > 1) {
        printf("hdr1\n");
        dump_hdr(hdr1);
    }
    if (verbose) printf("RUN test 1\n");

    // test
    redirected_stderr = redirect_stderr(tempfname);
    bool result_1 = (!sam_hdr_remove_except(hdr1, "RG", "ID", id_to_keep_1) &&
                     !sam_hdr_add_pg(hdr1, "samtools", "VN", samtools_version(),
                                     arg_list ? "CL": NULL,
                                     arg_list ? arg_list : NULL,
                                     NULL));
    flush_and_restore_stderr(orig_stderr, redirected_stderr);

    if (verbose) printf("END RUN test 1\n");
    if (verbose > 1) {
        printf("hdr1\n");
        dump_hdr(hdr1);
    }

    // check result
    res.l = 0;
    check = fopen(tempfname, "r");
    if ( result_1
        && check_test_1(hdr1)
        && kgetline(&res, (kgets_func *)fgets, check) < 0
        && (feof(check) || res.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) printf("FAIL test 1\n");
    }
    fclose(check);

    // teardown
    sam_hdr_destroy(hdr1);
    if (verbose) printf("END test 1\n");

    if (verbose) printf("BEGIN test 2\n");  // test eliminating a tag that is there
    sam_hdr_t* hdr2;
    const char* id_to_keep_2 = "fish";
    setup_test_2(&hdr2);
    if (verbose > 1) {
        printf("hdr2\n");
        dump_hdr(hdr2);
    }
    if (verbose) printf("RUN test 2\n");

    // test
    redirected_stderr = redirect_stderr(tempfname);
    bool result_2 = (!sam_hdr_remove_except(hdr2, "RG", "ID", id_to_keep_2) &&
            !sam_hdr_add_pg(hdr2, "samtools", "VN", samtools_version(),
                                    arg_list ? "CL": NULL,
                                    arg_list ? arg_list : NULL,
                                    NULL));
    flush_and_restore_stderr(orig_stderr, redirected_stderr);

    if (verbose) printf("END RUN test 2\n");
    if (verbose > 1) {
        printf("hdr2\n");
        dump_hdr(hdr2);
    }

    // check result
    res.l = 0;
    check = fopen(tempfname, "r");
    if ( result_2
        && check_test_2(hdr2)
        && kgetline(&res, (kgets_func *)fgets, check) < 0
        && (feof(check) || res.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) printf("FAIL test 2\n");
    }
    fclose(check);

    // teardown
    sam_hdr_destroy(hdr2);
    if (verbose) printf("END test 2\n");

    if (verbose) printf("BEGIN test 3\n");  // test eliminating a tag that is there
    sam_hdr_t* hdr3;
    setup_test_3(&hdr3);
    if (verbose > 1) {
        printf("hdr3\n");
        dump_hdr(hdr3);
    }
    if (verbose) printf("RUN test 3\n");

    // test
    redirected_stderr = redirect_stderr(tempfname);
    bool result_3 = (!sam_hdr_remove_except(hdr3, "RG", NULL, NULL) &&
            !sam_hdr_add_pg(hdr3, "samtools", "VN", samtools_version(),
                                    arg_list ? "CL": NULL,
                                    arg_list ? arg_list : NULL,
                                    NULL));
    flush_and_restore_stderr(orig_stderr, redirected_stderr);

    if (verbose) printf("END RUN test 3\n");
    if (verbose > 1) {
        printf("hdr3\n");
        dump_hdr(hdr3);
    }

    // check result
    res.l = 0;
    check = fopen(tempfname, "r");
    if ( result_3
        && check_test_3(hdr3)
        && kgetline(&res, (kgets_func *)fgets, check) < 0
        && (feof(check) || res.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) printf("FAIL test 3\n");
    }
    fclose(check);

    // teardown
    sam_hdr_destroy(hdr3);
    if (verbose) printf("END test 3\n");

    // Cleanup
    free(res.s);
    free(arg_list);
    remove(tempfname);
    if (failure > 0)
        fprintf(stderr, "%d failures %d successes\n", failure, success);
    close(orig_stderr);

    return (success == NUM_TESTS)? EXIT_SUCCESS : EXIT_FAILURE;
}
