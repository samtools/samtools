/*  test/split/test_filter_header_rg.c -- split test cases.

    Copyright (C) 2014 Genome Research Ltd.

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

#include "../../bam_split.c"
#include "../test.h"
#include <unistd.h>

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
    if (!lines2) return false;

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
    return res;
}

void setup_test_1(bam_hdr_t** hdr_in)
{
    *hdr_in = bam_hdr_init();
    const char *test1 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@RG\tID:fish\n";
    (*hdr_in)->text = strdup(test1);
    (*hdr_in)->l_text = strlen(test1);
}

bool check_test_1(const bam_hdr_t* hdr) {
    const char *test1_res =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@PG\tID:samtools\tPN:samtools\tVN:x.y.test\tCL:test_filter_header_rg foo bar baz\n";

    if (hdrcmp(hdr->text, test1_res)) {
        return false;
    }
    return true;
}

void setup_test_2(bam_hdr_t** hdr_in)
{
    *hdr_in = bam_hdr_init();
    const char *test2 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@RG\tID:fish\n";
    (*hdr_in)->text = strdup(test2);
    (*hdr_in)->l_text = strlen(test2);
}

bool check_test_2(const bam_hdr_t* hdr) {
    const char *test2_res =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@RG\tID:fish\n"
    "@PG\tID:samtools\tPN:samtools\tVN:x.y.test\tCL:test_filter_header_rg foo bar baz\n";

    if (hdrcmp(hdr->text, test2_res)) {
        return false;
    }
    return true;
}

int main(int argc, char *argv[])
{
    // test state
    const int NUM_TESTS = 2;
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
    bam_hdr_t* hdr1;
    const char* id_to_keep_1 = "1#2.3";
    setup_test_1(&hdr1);
    if (verbose > 1) {
        printf("hdr1\n");
        dump_hdr(hdr1);
    }
    if (verbose) printf("RUN test 1\n");

    // test
    redirected_stderr = redirect_stderr(tempfname);
    bool result_1 = filter_header_rg(hdr1, id_to_keep_1, arg_list);
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
    bam_hdr_destroy(hdr1);
    if (verbose) printf("END test 1\n");

    if (verbose) printf("BEGIN test 2\n");  // test eliminating a tag that is there
    bam_hdr_t* hdr2;
    const char* id_to_keep_2 = "fish";
    setup_test_2(&hdr2);
    if (verbose > 1) {
        printf("hdr2\n");
        dump_hdr(hdr2);
    }
    if (verbose) printf("RUN test 2\n");

    // test
    redirected_stderr = redirect_stderr(tempfname);
    bool result_2 = filter_header_rg(hdr2, id_to_keep_2, arg_list);
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
    bam_hdr_destroy(hdr2);
    if (verbose) printf("END test 2\n");


    // Cleanup
    free(res.s);
    free(arg_list);
    remove(tempfname);
    if (failure > 0)
        fprintf(stderr, "%d failures %d successes\n", failure, success);
    close(orig_stderr);

    return (success == NUM_TESTS)? EXIT_SUCCESS : EXIT_FAILURE;
}
