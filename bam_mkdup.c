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

/*
 * loop through bam
 *  collect a pair
 *  try insert pair into hash
 *   if it's already in there it's a duplicate so mark it if not insert it
 */

typedef struct possig_part {
	bool orient:1;
	uint32_t tid:31;
	uint32_t pos;
} possig_part_t;

typedef struct possig {
	union {
		struct {
			possig_part_t first;
			possig_part_t second;
		} field;
		__uint128_t bits;
	};
} possig_t;


// FIXME: write an alternative for platforms that don't have a 128 bit intrinsic
// FIXME: this is the 64bit hash find a good one for 128 bit numbers
#define possig_hash_func(key) (khint32_t)((key.bits)>>33^(key.bits)^(key.bits)<<11)

#define possig_hash_equal(a, b) ((a.bits) == (b.bits))

KHASH_INIT(sig, possig_t, char, 1, possig_hash_func, possig_hash_equal)

typedef struct parsed_opts {
	char* input_name;
	char* output_name;
} parsed_opts_t;

typedef struct state {
	samFile* fin;
	samFile* fout;
	bam_hdr_t* hin;
	bam_hdr_t* hout;
	kh_sig_t* hash; // hash of seen positions could replace with bloom filter?
} state_t;

static inline bool bam_to_possig_part(bam1_t* record, possig_part_t* part)
{
	if (record->core.tid > INT32_MAX) return false;
	part->orient = ((record->core.flag&BAM_FREVERSE) == BAM_FREVERSE);
	part->tid = record->core.tid;
	return true;
}

static inline bool bam_to_possig(bam1_t* first_pos, bam1_t* second_pos, possig_t* sig)
{
	// Get 5 prime coord of each part of reads
	int32_t fpos = (first_pos->core.flag&BAM_FREVERSE) ? bam_endpos(first_pos) : first_pos->core.pos;
	int32_t spos = (second_pos->core.flag&BAM_FREVERSE) ? bam_endpos(second_pos) : second_pos->core.pos;

	// Make the position signature with the first coordinate read in the first position
	if (first_pos->core.tid == second_pos->core.tid) {
		if (fpos < spos) {
			sig->field.first.pos = fpos;
			sig->field.second.pos = spos;
			return (bam_to_possig_part(first_pos, &sig->field.first) && bam_to_possig_part(second_pos, &sig->field.second));
		} else {
			sig->field.first.pos = spos;
			sig->field.second.pos = fpos;
			return (bam_to_possig_part(second_pos, &sig->field.first) && bam_to_possig_part(first_pos, &sig->field.second));
		}
	} else {
		if (first_pos->core.tid < second_pos->core.tid) {
			sig->field.first.pos = fpos;
			sig->field.second.pos = spos;
			return (bam_to_possig_part(first_pos, &sig->field.first) && bam_to_possig_part(second_pos, &sig->field.second));
		} else {
			sig->field.first.pos = spos;
			sig->field.second.pos = fpos;
			return (bam_to_possig_part(second_pos, &sig->field.first) && bam_to_possig_part(first_pos, &sig->field.second));
		}
	}
}

static bool loop(const state_t* state, bam1_t* first, bam1_t* second)
{
	// bam to position sig
	possig_t sig;
	if (!bam_to_possig(first, second, &sig)) { fprintf(stderr, "trace: unlikely tid\n"); return false; }
	
	int ret = 0;
	khiter_t val = kh_put(sig, state->hash, sig, &ret);
	if (ret == 0) {
		// Already in there ergo I am a duplicate
		first->core.flag |= BAM_FDUP;
		second->core.flag |= BAM_FDUP;
	} else {
		kh_value(state->hash, val) = 1;
	}
	return true;
}

/*
 * read a read
 *
 * if already have a read in prev
 *   if mismatch
 *      write prev
 *      if read is unpaired, unmapped, has unmapped mate, is secondary, qcfail or supplimentary
 *        write read
 *      else
 * 	   prev = read;
 *      endif
 *   elseif match
 *      process
 *      prev = NULL;
 *   end
 * else
 *   if read is unpaired, unmapped, has unmapped mate, is secondary, qcfail or supplimentary
 *     write read
 *   else
 *     put read in prev
 *   endif
 * endif
 */

static inline bool is_unprocessable( bam1_t* read )
{
	// if read is unpaired, unmapped, has unmapped mate, is secondary, qcfail or supplimentary
	return ((read->core.flag&(BAM_FUNMAP|BAM_FMUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FSUPPLEMENTARY)) != 0
	|| (read->core.flag&BAM_FPAIRED) != BAM_FPAIRED
	|| read->core.tid == -1
	|| read->core.pos == -1);
}

static bool process(const state_t* state)
{
	bam1_t* first = bam_init1();
	bam1_t* second = bam_init1();
	bool success = true;
	bool have_2nd = false;
	while (sam_read1(state->fin, state->hin, first) >= 0) {
		first->core.flag &= ~BAM_FDUP;  // Strip any existing duplicate flags
		bool unprocessable = is_unprocessable(first);
		if (have_2nd) {
			if (strcmp(bam_get_qname(first),bam_get_qname(second))) { // We've moved on
				sam_write1(state->fout, state->hout, second);
				if (unprocessable){
					sam_write1(state->fout, state->hout, first);
					have_2nd = false;
				} else {
					// put read in second
					bam1_t* tmp = second;
					second = first;
					first = tmp;
					// have_2nd remains true;
				}
			} else {
				if (unprocessable) {
					// this disrupts order but meh
					// FIXME: should implement stack to deal with this
					sam_write1(state->fout, state->hout, first);
				} else {
					if (!loop(state, first, second)) success = false;
					sam_write1(state->fout, state->hout, first);
					sam_write1(state->fout, state->hout, second);
					have_2nd = false;
				}
			}
		} else {
			if (unprocessable) {
				sam_write1(state->fout, state->hout, first);
			} else {
				// put read in second
				bam1_t* tmp = second;
				second = first;
				first = tmp;
				
				have_2nd = true;
			}
		}
	}
	if (have_2nd) {
		sam_write1(state->fout, state->hout, second);
	}

	bam_destroy1(first);
	bam_destroy1(second);
	return success;
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
	state->fin = sam_open(opts->input_name, "rb");
	state->hin = sam_hdr_read(state->fin);
	
	state->fout = sam_open(opts->output_name ? opts->output_name : "-", "wb");
	state->hout = bam_hdr_dup(state->hin);
	sam_hdr_write(state->fout, state->hout);
	
	state->hash = kh_init(sig);

	return true;
}

static void cleanup_state(state_t* state)
{
	sam_close(state->fout);
	sam_close(state->fin);
	kh_destroy(sig, state->hash);
}

int main_mkdup(int argc, char** argv)
{
	parsed_opts_t opts = {NULL, NULL};
	
	if ( !parse_args(argc, argv, &opts) ) return 1;
	
	int ret = 1;
	state_t state = {NULL, NULL, NULL, NULL, NULL};
	
	if ( init_state(&opts, &state)) {
		if (process(&state)) ret = 0;
	}
	cleanup_state(&state);

	return ret;
}
