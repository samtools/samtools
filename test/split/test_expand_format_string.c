/*  test/split/test_expand_format_string.c -- split format string test cases.

    Copyright (C) 2014-2015 Genome Research Ltd.

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
#include <stdlib.h>
#include <unistd.h>

void setup_test_1(sam_hdr_t** hdr_in)
{
    *hdr_in = sam_hdr_init();
    const char *test1 =
    "@HD\tVN:1.4\n"
    "@SQ\tSN:blah\n"
    "@RG\tID:fish\n";
    sam_hdr_add_lines(*hdr_in, test1, 0);
}

int main(int argc, char**argv)
{
    // test state
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
                printf(
                       "usage: test_expand_format_string [-v]\n\n"
                       " -v verbose output\n"
                       );
                break;
        }
    }


    // Setup stderr redirect
    kstring_t res = { 0, 0, NULL };
    FILE* orig_stderr = fdopen(dup(STDERR_FILENO), "a"); // Save stderr
    char* tempfname = (optind < argc)? argv[optind] : "test_expand_format_string.tmp";
    FILE* check = NULL;

    // setup
    if (verbose) printf("BEGIN test 1\n");  // default format string test
    const char* format_string_1 = "%*_%#.bam";
    const char* basename_1 = "basename";
    const char* rg_id_1 = "1#2.3";
    const int rg_idx_1 = 4;
    if (verbose > 1) {
        printf("format_string:%s\n"
               "basename:%s\n"
               "rg_id:%s\n"
               "rg_idx:%d\n", format_string_1, basename_1, rg_id_1, rg_idx_1);
    }
    if (verbose) printf("RUN test 1\n");

    // test
    xfreopen(tempfname, "w", stderr); // Redirect stderr to pipe
    char* output_1 = expand_format_string(format_string_1, basename_1, rg_id_1, rg_idx_1, NULL);
    fclose(stderr);

    if (verbose) printf("END RUN test 1\n");
    if (verbose > 1) {
        printf("format_string:%s\n"
               "basename:%s\n"
               "rg_id:%s\n"
               "rg_idx:%d\n", format_string_1, basename_1, rg_id_1, rg_idx_1);
    }

    // check result
    res.l = 0;
    check = fopen(tempfname, "r");
    if (output_1 != NULL && !strcmp(output_1, "basename_4.bam")
        && kgetline(&res, (kgets_func *)fgets, check) < 0
        && (feof(check) || res.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) printf("FAIL test 1\n");
    }
    fclose(check);

    // teardown
    free(output_1);
    if (verbose) printf("END test 1\n");

    // Cleanup test harness
    free(res.s);
    remove(tempfname);
    if (failure > 0)
        fprintf(orig_stderr, "%d failures %d successes\n", failure, success);
    fclose(orig_stderr);

    return (success == NUM_TESTS)? EXIT_SUCCESS : EXIT_FAILURE;
}
