/*  test/merge/test_pretty_header.c -- header test harness.

    Copyright (C) 2013 Genome Research Ltd.

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

void setup_test_1(char** input) {
    *input = strdup(
    "@HD\n"
    "@SQ\n"
    "@RG\n"
    "@PG\n"
    "@CO\n");
}

bool check_test_1(char* input) {
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

    if (verbose) printf("BEGIN test 1\n");
    // setup
    char* input;
    setup_test_1(&input);
    // test
    if (verbose > 1) {
        printf("input:\n%s",input);
    }
    if (verbose) printf("RUN test 1\n");
    pretty_header(&input, strlen(input));
    if (verbose) printf("END RUN test 1\n");
    if (verbose > 1) {
        printf("input:\n%s",input);
    }
    if (check_test_1(input)) { ++success; } else { ++failure; }
    // teardown
    free(input);
    if (verbose) printf("END test 1\n");

    if (success == NUM_TESTS) {
        return 0;
    } else {
        fprintf(stderr, "%d failures %d successes\n", failure, success);
        return 1;
    }
}
