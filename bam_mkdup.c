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


#include <htslib/sam.h>
#include <htslib/khash.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdio.h>
#include "bam_mkdup.h"

static const int SO_UNKNOWN = -1;
static const int SO_COORDINATE = 0;
static const int SO_NAME = 1;

typedef struct parsed_opts {
	char* input_name;
	char* output_name;
} parsed_opts_t;

static int detect_sort_order(const bam_hdr_t* test)
{
	char* first_n_ptr = strchr(test->text, '\n');
	if ( test->l_text < 4 || strncmp(test->text,"@HD",3) || first_n_ptr == NULL) return SO_UNKNOWN;
	
	int first_n = first_n_ptr - test->text;
	char* parse = strndup(test->text, first_n);
	int ret = SO_UNKNOWN;
	
	if (strstr(parse, "\tSO:coordinate")) { ret = SO_COORDINATE; }
	else if (strstr(parse, "\tSO:queryname")) { ret = SO_NAME; }
	free(parse);

	return ret;
}

static void usage(FILE* where)
{
	fprintf(where, "Usage: samtools mkdup <input.bam> [<output.bam>]\n\n"
			"Marks duplicates within a file.\n"
			"If - is used for input.bam stdin will be used.\n"
			"If output.bam is not specified stdout will be assumed.\n"
			);
}

static bool parse_args( int argc, char** argv, parsed_opts_t* opts )
{
	// Check number of input arguments, minimum 1, maximum 2
	if (argc == 1) { usage(stdout); return false; }
	if (argc < 2 || argc > 3) { usage(stderr); return false; }
	
	opts->input_name = argv[1];
	if ( argc == 3 ) opts->output_name = argv[2];
	
	return true;
}

static bool init_state(const parsed_opts_t* opts, state_t* state)
{
	state->fin = sam_open(opts->input_name, "r");
	state->hin = sam_hdr_read(state->fin);
	
	if ((state->so = detect_sort_order(state->hin)) == SO_UNKNOWN) {
		fprintf(stderr, "[init_state] Sort order unknown, cannot proceed.\n");
		return false;
	}
	
	state->fout = sam_open(opts->output_name ? opts->output_name : "-", "wb");
	state->hout = bam_hdr_dup(state->hin);
	if (sam_hdr_write(state->fout, state->hout) < 0) {
		fprintf(stderr, "[init_state] Unable to write header.\n");
		return false;
	}
	
	if ( state->so == SO_NAME ) state->hash = kh_init(sig);

	return true;
}

static void cleanup_state(state_t* state)
{
	if (state) {
		sam_close(state->fout);
		sam_close(state->fin);
		kh_destroy(sig, state->hash);
	}
}

int main_mkdup(int argc, char** argv)
{
	parsed_opts_t opts = {NULL, NULL};
	
	if ( !parse_args(argc, argv, &opts) ) return 1;
	
	int ret = 1;
	state_t state = {NULL, NULL, NULL, NULL, 0, NULL};
	
	if ( init_state(&opts, &state)) {
		if (state.so == SO_NAME) {
			if (process_namesorted(&state)) ret = 0;
		} else {
			if (process_coordsorted(&state, opts.input_name)) ret = 0;
		}
	}
	cleanup_state(&state);

	return ret;
}
