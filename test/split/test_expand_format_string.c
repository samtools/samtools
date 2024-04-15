/*  test/split/test_expand_format_string.c -- split format string test cases.

    Copyright (C) 2014-2015,2024 Genome Research Ltd.

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

typedef struct {
    char *tempfname;
    kstring_t res;
    int verbose;
    int test_num;
    int success;
    int failure;
} TestState;

void run_test(TestState *state, const char* format_string, const char* basename,
              const char* rg_id, const int rg_idx, const int padding,
              const htsFormat *format, const char *expected) {
    FILE* check = NULL;

    ++state->test_num;
    if (state->verbose) printf("BEGIN test %d\n", state->test_num);
    if (state->verbose > 1) {
        printf("format_string:%s\n"
               "basename:%s\n"
               "rg_id:%s\n"
               "rg_idx:%d\n", format_string, basename, rg_id, rg_idx);
    }
    if (state->verbose) printf("RUN test %d\n", state->test_num);

    // test
    xfreopen(state->tempfname, "w", stderr); // Redirect stderr to pipe
    char* output = expand_format_string(format_string, basename,
                                        rg_id, rg_idx, padding, format);
    fclose(stderr);

    if (state->verbose) printf("END RUN test %d\n", state->test_num);
    if (state->verbose > 1) {
        printf("got: \"%s\" expected: \"%s\"\n", output ? output : "(null)",
               expected ? expected : "(null)");
    }
    if (state->verbose > 2) {
        printf("format_string:%s\n"
               "basename:%s\n"
               "rg_id:%s\n"
               "rg_idx:%d\n", format_string, basename, rg_id, rg_idx);
    }

    // check result
    state->res.l = 0;
    check = fopen(state->tempfname, "r");
    if (expected != NULL) {
        // Good input, should produce a file name
        if (check && output != NULL && !strcmp(output, expected)
            && kgetline(&state->res, (kgets_func *)fgets, check) < 0
            && (feof(check) || state->res.l == 0)) {
            ++state->success;
        } else {
            ++state->failure;
            if (state->verbose) printf("FAIL test %d\n", state->test_num);
        }
    } else {
        // Bad input, should return NULL and print an error message
        if (check && output == NULL
            && kgetline(&state->res, (kgets_func *)fgets, check) == 0
            && state->res.l > 0) {
            ++state->success;
        } else {
            ++state->failure;
            if (state->verbose) printf("FAIL test %d\n", state->test_num);
        }
    }
    fclose(check);

    // teardown
    free(output);
    if (state->verbose) printf("END test %d\n", state->test_num);
    return;
}

int main(int argc, char**argv)
{
    // test state
    TestState state = { NULL, KS_INITIALIZE, 0, 0, 0, 0};
    const static htsFormat sam_fmt  = { sequence_data, sam  };
    const static htsFormat bam_fmt  = { sequence_data, bam  };
    const static htsFormat cram_fmt = { sequence_data, cram };

    int getopt_char;
    while ((getopt_char = getopt(argc, argv, "v")) != -1) {
        switch (getopt_char) {
            case 'v':
                ++state.verbose;
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
    FILE* orig_stderr = fdopen(dup(STDERR_FILENO), "a"); // Save stderr
    state.tempfname = (optind < argc)? argv[optind] : "test_expand_format_string.tmp";

    // default format string test
    run_test(&state, "%*_%#.%.", "basename", "1#2.3", 4, 0, &bam_fmt,
             "basename_4.bam");

    // default for non-RG tags
    run_test(&state, "%*_%!.%.", "basename", "1#2.3", 4, 0, &bam_fmt,
             "basename_1#2.3.bam");

    // alternate file types
    run_test(&state, "%*_%#.%.", "basename", "1#2.3", 4, 0, &sam_fmt,
             "basename_4.sam");
    run_test(&state, "%*_%#.%.", "basename", "1#2.3", 4, 0, &cram_fmt,
             "basename_4.cram");

    // zero padding
    run_test(&state, "%*_%#.%.", "basename", "1#2.3", 4, 5, &bam_fmt,
             "basename_00004.bam");

    // percent sign
    run_test(&state, "%%%*_%#.%.%%", "basename", "1#2.3", 4, 0, &bam_fmt,
             "%basename_4.bam%");

    // Bad input - single percent sign at end
    run_test(&state, "%%%*_%#.%.%", "basename", "1#2.3", 4, 0, &bam_fmt,
             NULL);

    // Bad input - invalid format character
    run_test(&state, "%s_%#.%.", "basename", "1#2.3", 4, 0, &bam_fmt,
             NULL);


    // Cleanup test harness
    free(state.res.s);
    remove(state.tempfname);
    if (state.failure > 0)
        fprintf(orig_stderr, "%d failures %d successes\n",
                state.failure, state.success);
    fclose(orig_stderr);

    return (state.success == state.test_num)? EXIT_SUCCESS : EXIT_FAILURE;
}
