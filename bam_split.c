/* The MIT License
 *
 * Copyright (c) 2013, 2014 Genome Research Limited.
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

#define _GNU_SOURCE

#include <htslib/sam.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>
#include <regex.h>
#include <htslib/khash.h>

KHASH_MAP_INIT_STR(c2i, int)

struct parsed_opts {
	char* merged_output_name;
	char* unaccounted_header_name;
	char* unaccounted_name;
	char* output_prefix;
	bool verbose;
};

typedef struct parsed_opts parsed_opts_t;

struct state {
	samFile* merged_input_file;
	bam_hdr_t* merged_input_header;
	samFile* unaccounted_file;
	bam_hdr_t* unaccounted_header;
	size_t output_count;
	char** rg_id;
	samFile** rg_output_file;
	bam_hdr_t** rg_output_header;
	kh_c2i_t* rg_hash;
};

typedef struct state state_t;

static void cleanup_state(state_t* status);
static void cleanup_opts(parsed_opts_t* opts);

static void usage(bool error)
{
	FILE* write_to = error ? stderr : stdout;
	fprintf(write_to,
			"samtools split\n\n"
			"usage: samtools split [-v] <merged.bam> <unaccounted.bam>[:<unaccounted_header.sam>] <prefix>\n"
			"Where:\n"
			"<merged.bam> is your file to separate out into individual read groups."
			"<unaccounted.bam> is where to put reads without a read group tag.\n"
			"<:unaccounted_header.sam> is an optional header for <unaccounted.bam>.\n"
			"-v    verbose output.\n"
			);
}

// Takes the command line options and turns them into something we can understand
static parsed_opts_t* parse_args(int argc, char** argv)
{
	const char* optstring = "v";
	
	parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
	if (! retval ) return NULL;
	
	int opt;
	while ((opt = getopt(argc, argv, optstring)) != -1) {
		switch (opt) {
		case 'v':
			retval->verbose = true;
			break;
		default:
			usage(false);
			free(retval);
			return NULL;
		}
	}
	
	if (argc < 3) {
		usage(true);
		free(retval);
		return NULL;
	}

	retval->merged_output_name = strdup(argv[1]);
	char* sep = strdup(argv[2]);
	retval->unaccounted_name = strsep(&sep, ":");
	retval->unaccounted_header_name = sep;
	retval->output_prefix = strdup(argv[3]);

	return retval;
}

// Parse the header, count the number of RG tags and return a list of their names
static bool count_RG(bam_hdr_t* hdr, size_t* count, char*** output_name)
{
	if (hdr->l_text < 3 ) {
		*count = 0;
		*output_name = NULL;
		return true;
	}
	char* input = strndup(hdr->text, hdr->l_text);
	if ( input == NULL ) {
		return false;
	}
	
	//////////////////////////////////////////
	// First stage count number of @RG tags //
	//////////////////////////////////////////
	char* pointer = input;
	size_t n_rg = 0;
	// Guard against rare case where @RG is first header line
	// This shouldn't happen but could where @HD is omitted
	if (pointer[0] == '@' && pointer[1] == 'R' && pointer[2] == 'G' ) {
		++n_rg;
		pointer += 3;
	}
	char* line;
	while ((line = strstr(pointer, "\n@RG")) != NULL) {
		++n_rg;
		pointer = line + 1;
	}
		
	//////////////////////////////////
	// Second stage locate @RG ID's //
	//////////////////////////////////
	char** names = (char**)calloc(sizeof(char*), n_rg);
	size_t next = 0;
	
	regex_t rg_finder;
	if (regcomp(&rg_finder, "^@RG.*\tID:([!-)+-<>-~][ !-~]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE) != 0) {
		free(input);
		free(names);
		return false;
	}
	regmatch_t* matches = (regmatch_t*)calloc(sizeof(regmatch_t),2);
	int error;
	char* begin = input;

	while ((error = regexec(&rg_finder, begin, 2, matches, 0)) == 0) {
		names[next++] = strndup(begin+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so);
		begin += matches[0].rm_eo;
	}

	if (error != REG_NOMATCH) {
		// cleanup
		regfree(&rg_finder);
		free(matches);
		free(names);
		free(input);
		return false;
	}
	free(matches);

	// return results
	*count = n_rg;
	*output_name = names;
	regfree(&rg_finder);
	free(input);
	return true;
}

// Set the initial state
static state_t* init(parsed_opts_t* opts)
{
	state_t* retval = malloc(sizeof(state_t));
	if (!retval) {
		fprintf(stderr, "Out of memory");
		return NULL;
	}

	retval->merged_input_file = sam_open(opts->merged_output_name, "rb");
	if (!retval->merged_input_file) {
		fprintf(stderr, "Could not open header file (%s)", opts->merged_output_name);
		return NULL;
	}
	retval->merged_input_header = sam_hdr_read(retval->merged_input_file);
	if (opts->unaccounted_header_name) {
		samFile* hdr_load = sam_open(opts->unaccounted_header_name, "r");
		if (!hdr_load) {
			fprintf(stderr, "Could not open unaccounted header file (%s)", opts->unaccounted_header_name);
			return NULL;
		}
		retval->unaccounted_header = sam_hdr_read(hdr_load);
		sam_close(hdr_load);
	} else {
		retval->unaccounted_header = bam_hdr_dup(retval->merged_input_header);
	}

	retval->unaccounted_file = sam_open(opts->unaccounted_name, "wb");
	if (retval->unaccounted_file == NULL) {
		fprintf(stderr, "Could not open unaccounted output file: %s\r\n", opts->unaccounted_name);
		cleanup_state(retval);
		return NULL;
	}

	// Open output files for RGs
	if (!count_RG(retval->merged_input_header, &retval->output_count, &retval->rg_id)) return NULL;
	if (opts->verbose) fprintf(stderr, "@RG's found %zu\n",retval->output_count);

	retval->rg_output_file = (samFile**)calloc(retval->output_count, sizeof(samFile*));
	retval->rg_output_header = (bam_hdr_t**)calloc(retval->output_count, sizeof(bam_hdr_t*));
	retval->rg_hash = kh_init_c2i();
	if (!retval->rg_output_file || !retval->rg_output_header) {
		fprintf(stderr, "Could not allocate memory for output file array. Out of memory?");
		free(retval);
		return NULL;
	}
	size_t i;
	for (i = 0; i < retval->output_count; i++) {
		char* output_filename = NULL;
		if (asprintf( &output_filename, "%s%s.bam",opts->output_prefix,retval->rg_id[i]) == -1) {
			fprintf(stderr, "Could not allocate memory for output filename string.\r\n");
			return NULL;
		}
		retval->rg_output_file[i] = sam_open(output_filename, "wb");
		if (retval->rg_output_file[i] == NULL) {
			fprintf(stderr, "Could not open output file: %s\r\n", output_filename);
			return NULL;
		}

		// Record index in hash
		int ret;
		khiter_t iter = kh_put_c2i(retval->rg_hash, retval->rg_id[i], &ret);
		kh_val(retval->rg_hash,iter) = i;

		free(output_filename);
		retval->rg_output_header[i] = bam_hdr_dup(retval->merged_input_header);
	}

	return retval;
}

static bool split(state_t* state)
{
	if (sam_hdr_write(state->unaccounted_file, state->unaccounted_header) != 0) {
		fprintf(stderr, "Could not write output file header");
		return false;
	}
	size_t i;
	for (i = 0; i < state->output_count; i++) {
		if (sam_hdr_write(state->rg_output_file[i], state->rg_output_header[i]) != 0) {
			fprintf(stderr, "Could not write output file header");
			return false;
		}
	}

	bam1_t* file_read = bam_init1();
	// Read the first record
	if (sam_read1(state->merged_input_file, state->merged_input_header, file_read) < 0) {
		// Nothing more to read?  Ignore this file
		bam_destroy1(file_read);
		file_read = NULL;
	}

	while (file_read != NULL) {
		// Get RG tag from read and look it up in hash to find file to output it to
		uint8_t* tag = bam_aux_get(file_read, "RG");
		khiter_t iter;
		if ( tag != NULL ) {
			char* rg = bam_aux2Z(tag);
			iter = kh_get_c2i(state->rg_hash, rg);
		} else {
			iter = kh_end(state->rg_hash);
		}
		
		// Write the read out to correct file
		if (iter != kh_end(state->rg_hash)) {
			// if found write to the appropriate untangled bam
			int i = kh_val(state->rg_hash,iter);
			sam_write1(state->rg_output_file[i], state->rg_output_header[i], file_read);
		} else {
			// otherwise write to the unaccounted bam
			sam_write1(state->unaccounted_file, state->unaccounted_header, file_read);
		}

		// Replace written read with the next one to process
		if (sam_read1(state->merged_input_file, state->merged_input_header, file_read) < 0) {
			// Nothing more to read?  Ignore this file in future
			bam_destroy1(file_read);
			file_read = NULL;
		}
	}

	return true;
}

static void cleanup_state(state_t* status)
{
	sam_close(status->unaccounted_file);
	sam_close(status->merged_input_file);
	size_t i;
	for (i = 0; i < status->output_count; i++) {
		bam_hdr_destroy(status->rg_output_header[i]);
		sam_close(status->rg_output_file[i]);
	}
	bam_hdr_destroy(status->merged_input_header);
	free(status->rg_output_header);
	free(status->rg_output_file);
	kh_destroy_c2i(status->rg_hash);
}

static void cleanup_opts(parsed_opts_t* opts)
{
	free(opts->merged_output_name);
	free(opts->unaccounted_header_name);
	free(opts->unaccounted_name);
	free(opts->output_prefix);
}

int main_split(int argc, char** argv)
{
	parsed_opts_t* opts = parse_args(argc, argv);
	if (!opts ) return 1;
	state_t* status = init(opts);
	if (!status) return 1;

	if (!split(status)) return 1;

	cleanup_state(status);
	cleanup_opts(opts);
	
	return 0;
}
