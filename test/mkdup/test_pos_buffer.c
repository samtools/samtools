/*  test/merge/test_pos_buff.c -- position buffer test harness.

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

#include "../../pos_buffer.c"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

void dump_rv(const read_vector_t rv) {
    printf("orient = %u, tid = %u, pos = %u\n",rv.orient, rv.tid, rv.pos);
}

void setup_test_1(read_vector_t* rv) {
    const read_vector_t rv_in[] = {
        { .orient = false, .tid = 0, .pos = 0 },
        { .orient = true, .tid = 0, .pos = 1 }
    };
    memcpy(rv, rv_in, sizeof(rv_in));
}

int main(int argc, char**argv)
{
    // test state
    const int NUM_TESTS = 2;
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

    // Setup stderr redirect
    size_t len = 0;
    char* res = NULL;
    FILE* orig_stderr = fdopen(dup(STDERR_FILENO), "a"); // Save stderr
    char* tempfname = (optind < argc)? argv[optind] : "test_pos_buffer.tmp";
    FILE* check = NULL;

    // setup
    if (verbose) printf("BEGIN test 1\n");  // Insert test
    read_vector_t rv1[2];
    pos_buffer_t* buf1 = pos_buffer_init();
    setup_test_1(rv1);
    if (verbose > 1) {
        printf("rv[0]:");
        dump_rv(rv1[0]);
        printf("rv[1]:");
        dump_rv(rv1[1]);
    }
    if (verbose) printf("RUN test 1\n");

    // test
    xfreopen(tempfname, "w", stderr); // Redirect stderr to pipe
    char* result = pos_buffer_insert(buf1, rv1[0], rv1[1], 20, "r1", 0, 1);
    fclose(stderr);

    if (verbose) printf("END RUN test 1\n");
    if (verbose > 1) {
        printf("rv[0]:");
        dump_rv(rv1[0]);
        printf("rv[1]:");
        dump_rv(rv1[1]);
    }

    // check result
    check = fopen(tempfname, "r");
    if ( (getline(&res, &len, check) == -1 ) &&
        (feof(check) || (res && !strcmp("",res))) &&
        result == NULL) {
        ++success;
    } else {
        ++failure;
        if (verbose) printf("FAIL test 1\n");
    }
    fclose(check);

    // teardown
    pos_buffer_destroy(buf1);
    if (verbose) printf("END test 1\n");

    // setup
    if (verbose) printf("BEGIN test 2\n");  // Insert test
    read_vector_t rv2[2];
    pos_buffer_t* buf2 = pos_buffer_init();
    setup_test_1(rv2);
    if (verbose > 1) {
        printf("rv[0]:");
        dump_rv(rv2[0]);
        printf("rv[1]:");
        dump_rv(rv2[1]);
    }
    if (verbose) printf("RUN test 2\n");
    
    // test
    xfreopen(tempfname, "w", stderr); // Redirect stderr to pipe
    char* result_1 = pos_buffer_insert(buf2, rv2[0], rv2[1], 20, "r1", 0, 1);
    char* result_2 = pos_buffer_insert(buf2, rv2[0], rv2[1], 19, "r2", 0, 1);
    fclose(stderr);
    
    if (verbose) printf("END RUN test 2\n");
    if (verbose > 1) {
        printf("rv[0]:");
        dump_rv(rv2[0]);
        printf("rv[1]:");
        dump_rv(rv2[1]);
        printf("result_1:%s\nresult_2:%s", result_1, result_2);
    }
    
    // check result
    check = fopen(tempfname, "r");
    if ( (getline(&res, &len, check) == -1 ) &&
        (feof(check) || (res && !strcmp("",res))) &&
        result_1 == NULL &&
        result_2 != NULL) {
        ++success;
    } else {
        ++failure;
        if (verbose) printf("FAIL test 2\n");
    }
    fclose(check);
    
    // teardown
    pos_buffer_destroy(buf2);
    if (verbose) printf("END test 2\n");
    
    // Cleanup
    free(res);
    remove(tempfname);
    if (failure > 0)
        fprintf(orig_stderr, "%d failures %d successes\n", failure, success);
    fclose(orig_stderr);

    return (success == NUM_TESTS)? EXIT_SUCCESS : EXIT_FAILURE;
}
