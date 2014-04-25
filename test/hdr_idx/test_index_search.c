/* The MIT License
 *
 * Copyright (c) 2014 Genome Research Limited.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Contact: Martin Pollard <mp15@sanger.ac.uk> */
/* Created by Martin Pollard on 25/04/2014. */

#include "../../hdr_idx_priv.c"
#include "../test.h"
#include <unistd.h>
#include <assert.h>

void setup_test_1(bam_hdr_t** hdr_in, bam1_t** read_in)
{
	*hdr_in = bam_hdr_init();
	const char *test1 =
	"@HD\tVN:1.4\n"
	"@SQ\tSN:blah\n"
	"@RG\tID:fish\tLB:1234\n";
	(*hdr_in)->text = strdup(test1);
	(*hdr_in)->l_text = strlen(test1);
	*(read_in) = bam_init1();
	bam_aux_append(*(read_in), "RG", 'Z', 5, (uint8_t*)"fish");
}

void setup_test_2(bam_hdr_t** hdr_in, bam1_t** read_in)
{
	*hdr_in = bam_hdr_init();
	const char *test1 =
	"@HD\tVN:1.4\n"
	"@SQ\tSN:blah\n"
	"@RG\tID:fish\tLB:1234\n";
	(*hdr_in)->text = strdup(test1);
	(*hdr_in)->l_text = strlen(test1);
	*(read_in) = bam_init1();
}

void setup_test_3(bam_hdr_t** hdr_in, bam1_t** read_in)
{
	*hdr_in = bam_hdr_init();
	const char *test1 =
	"@HD\tVN:1.4\n"
	"@SQ\tSN:blah\n"
	"@RG\tLB:1234\tID:fish\n";
	(*hdr_in)->text = strdup(test1);
	(*hdr_in)->l_text = strlen(test1);
	*(read_in) = bam_init1();
	bam_aux_append(*(read_in), "RG", 'Z', 5, (uint8_t*)"fish");
}

int main(int argc, char**argv)
{
	// test state
	const int NUM_TESTS = 3;
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
					   "usage: test_index_search [-v]\n\n"
					   " -v verbose output\n"
					   );
				break;
		}
	}

	// Setup stderr redirect
	size_t len = 0;
	char* res = NULL;
	FILE* orig_stderr = fdopen(dup(STDERR_FILENO), "a"); // Save stderr
	char* tempfname = (optind < argc)? argv[optind] : "test_count_rg.tmp";
	FILE* check = NULL;

	// setup
	if (verbose) printf("BEGIN test 1\n");  // test init of library index under normal conditions
	bam1_t* search_1;
	bam_hdr_t* hdr1;
	setup_test_1(&hdr1, &search_1);
	library_index_t* index_1 = bam_library_index_init(hdr1);
	assert(index_1);
	if (verbose > 1) {
		printf("hdr1\n");
		dump_hdr(hdr1);
		printf("search_1\n");
		dump_read(search_1);
	}
	if (verbose) printf("RUN test 1\n");

	xfreopen(tempfname, "w", stderr); // Redirect stderr to pipe
	const char* result_1 = bam_search_library_index(index_1, search_1);
	fclose(stderr);

	if (verbose) printf("END RUN test 1\n");
	if (verbose > 1) {
		printf("hdr1\n");
		dump_hdr(hdr1);
		printf("result_1: \"%s\"\n", result_1);
		printf("search_1\n");
		dump_read(search_1);
	}

	// check result
	check = fopen(tempfname, "r");
	if (result_1 && !strcmp(result_1, "1234")
		&& (getline(&res, &len, check) == -1)
		&& (feof(check) || (res && !strcmp("",res)))) {
		++success;
	} else {
		++failure;
		if (verbose) printf("FAIL test 1\n");
	}
	fclose(check);

	// teardown
	bam_library_index_destroy(index_1);
	bam_hdr_destroy(hdr1);
	bam_destroy1(search_1);
	if (verbose) printf("END test 1\n");

	// setup
	if (verbose) printf("BEGIN test 2\n");  // test of read with no RG tag
	bam1_t* search_2;
	bam_hdr_t* hdr2;
	setup_test_2(&hdr2, &search_2);
	library_index_t* index_2 = bam_library_index_init(hdr2);
	if (verbose > 1) {
		printf("hdr2\n");
		dump_hdr(hdr2);
		printf("search_2\n");
		dump_read(search_2);
	}
	if (verbose) printf("RUN test 2\n");

	// test init of library normal conditions
	xfreopen(tempfname, "w", stderr); // Redirect stderr to pipe
	const char* result_2 = bam_search_library_index(index_2, search_2);
	fclose(stderr);

	if (verbose) printf("END RUN test 2\n");
	if (verbose > 1) {
		printf("hdr2\n");
		dump_hdr(hdr2);
		printf("result_2: \"%s\"\n", result_2);
		printf("search_2\n");
		dump_read(search_2);
	}

	// check result
	check = fopen(tempfname, "r");
	if (!result_2
		&& (getline(&res, &len, check) == -1)
		&& (feof(check) || (res && !strcmp("",res)))) {
		++success;
	} else {
		++failure;
		if (verbose) printf("FAIL test 2\n");
	}
	fclose(check);

	// teardown
	bam_library_index_destroy(index_2);
	bam_hdr_destroy(hdr2);
	bam_destroy1(search_2);
	if (verbose) printf("END test 2\n");

	// setup
	if (verbose) printf("BEGIN test 3\n");  // test search of library index under normal conditions with ID and LB in diff order
	bam1_t* search_3;
	bam_hdr_t* hdr3;
	setup_test_3(&hdr3, &search_3);
	library_index_t* index_3 = bam_library_index_init(hdr3);
	assert(index_3);
	if (verbose > 1) {
		printf("hdr3\n");
		dump_hdr(hdr3);
		printf("search_3\n");
		dump_read(search_3);
	}
	if (verbose) printf("RUN test 3\n");

	xfreopen(tempfname, "w", stderr); // Redirect stderr to pipe
	const char* result_3 = bam_search_library_index(index_3, search_3);
	fclose(stderr);

	if (verbose) printf("END RUN test 3\n");
	if (verbose > 1) {
		printf("hdr3\n");
		dump_hdr(hdr3);
		printf("result_3: \"%s\"\n", result_3);
		printf("search_3\n");
		dump_read(search_3);
		khint_t k;
		printf("index_3\n");
		for (k = kh_begin(index_3); k < kh_end(index_3); ++k) {
			if (kh_exist(index_3, k)) {
				printf("(%s)->(%s)\n", kh_key(index_3, k), kh_value(index_3, k));
			}
		}
	}

	// check result
	check = fopen(tempfname, "r");
	if (result_3 && !strcmp(result_3, "1234")
		&& (getline(&res, &len, check) == -1)
		&& (feof(check) || (res && !strcmp("",res)))) {
		++success;
	} else {
		++failure;
		if (verbose) printf("FAIL test 3\n");
	}
	fclose(check);

	// teardown
	bam_library_index_destroy(index_3);
	bam_hdr_destroy(hdr3);
	bam_destroy1(search_3);
	if (verbose) printf("END test 3\n");

	// Cleanup
	free(res);
	remove(tempfname);

	if (NUM_TESTS == success) {
		return 0;
	} else {
		fprintf(orig_stderr, "%d failures %d successes\n", failure, success);
		return 1;
	}
}
