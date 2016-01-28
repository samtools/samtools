/*  test/tview/test_get_rg_sample.c -- tview test cases.

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

#include "../../bam_tview.c"
#include <stdbool.h>

const char header_1[] =
"@HD    VN:1.4  SO:undefined\n"
"@SQ    SN:dummy\n"
"@RG    ID:blah SM:po\n";

void setup_test_1(char** header)
{
    *header = strdup(header_1);
}

khash_t(kh_rg)* run_test_1(char* header)
{
    khash_t(kh_rg)* test_result = get_rg_sample(header,"po");
    return test_result;
}

bool check_test_1(khash_t(kh_rg)* test_result, char* header)
{
    if (strcmp(header_1, header)) return false;
    // test blah is in there
    if (kh_get(kh_rg, test_result, "blah") == kh_end(test_result))
    {
        return false;
    }
    return true;
}

void teardown_1(khash_t(kh_rg)* test_result, char* header)
{
    free(header);
}

int main(int argc, char** argv)
{
    const int NUM_TESTS = 1;
    int success = 0;
    int failure = 0;

    char* test_header_1;
    setup_test_1(&test_header_1);
    khash_t(kh_rg)* test_result_1 = run_test_1(test_header_1);
    if (!check_test_1(test_result_1, test_header_1))
        failure++;
    else
        success++;
    teardown_1(test_result_1, test_header_1);

    if (success == NUM_TESTS) {
        return 0;
    } else {
        fprintf(stderr, "%d failures %d successes\n", failure, success);
        return 1;
    }
}
