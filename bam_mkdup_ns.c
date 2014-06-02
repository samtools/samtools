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


typedef struct {
	size_t n, max;
	bam1_t **a;
} tmp_stack_t;

// Adds a read to the stack
static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
	if (stack->n == stack->max) {
		stack->max = stack->max? stack->max<<1 : 0x10000;
		stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
	}
	stack->a[stack->n++] = b;
}
// Writes out all reads in the stack to file out and deallocates the memory they're using
static inline void flush_stack_to_disk(tmp_stack_t *stack, samFile* out, const bam_hdr_t *h, const bool is_dup)
{
	int i;
	for (i = 0; i < stack->n; ++i) {
		if (is_dup) {
			stack->a[i]->core.flag |= BAM_FDUP;
		} else {
			stack->a[i]->core.flag &= ~BAM_FDUP;
		}
		sam_write1(out, h, stack->a[i]);
		bam_destroy1(stack->a[i]);
	}
	stack->n = 0;
}


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

struct possig {
	union {
		struct {
			possig_part_t first;
			possig_part_t second;
		} field;
		__uint128_t bits;
	};
};


// FIXME: write an alternative for platforms that don't have a 128 bit intrinsic
// FIXME: this is the 64bit hash find a good one for 128 bit numbers
#define possig_hash_func(key) (khint32_t)((key.bits)>>33^(key.bits)^(key.bits)<<11)

#define possig_hash_equal(a, b) ((a.bits) == (b.bits))

__KHASH_IMPL(sig, , possig_t, char, 1, possig_hash_func, possig_hash_equal)

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

static bool loop(const state_t* state, bam1_t* first, bam1_t* second, bool* is_dup)
{
	// bam to position sig
	possig_t sig;
	if (!bam_to_possig(first, second, &sig)) { fprintf(stderr, "trace: unlikely tid\n"); return false; }
	
	int ret = 0;
	khiter_t val = kh_put(sig, state->hash, sig, &ret);
	if (ret == 0) {
		// Already in there ergo I am a duplicate
		*is_dup = true;
	} else {
		kh_value(state->hash, val) = 1;
		*is_dup = false;
	}
	return true;
}

/*
 * read a read
 *
 * TODO: rewrite for stack based implementation
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

bool process_namesorted(const state_t* state)
{
	tmp_stack_t stack = { 0, 0, NULL };
	bam1_t* read = bam_init1();
	bam1_t* second = NULL;
	bool success = true;
	
	while (sam_read1(state->fin, state->hin, read) >= 0) {
		bool unprocessable = is_unprocessable(read);
		if (second != NULL) {
			if (strcmp(bam_get_qname(read),bam_get_qname(second))) {
				// We've moved on but not found a match to our previous processable read
				// TODO: process as single read here
				flush_stack_to_disk(&stack, state->fout, state->hout, false);
				stack_insert(&stack, read);
				if (unprocessable){
					second = NULL;
				} else {
					// put read in second
					second = read;
					// have_2nd remains true;
				}
				read = bam_init1();
			} else {
				if (!unprocessable) {
					bool is_dup = false;
					// We have a viable pair, check if we have duplicate
					if (!loop(state, read, second, &is_dup)) success = false;
					stack_insert(&stack, read);
					read = bam_init1();
					int x;
					// short circuit reading the rest of reads of same name
					while ((x = sam_read1(state->fin, state->hin, read)) >= 0 && !strcmp(bam_get_qname(read),bam_get_qname(second)) ) {
						stack_insert(&stack, read);
						read = bam_init1();
					}
					second = NULL;
					flush_stack_to_disk(&stack, state->fout, state->hout, is_dup);
					if (x >= 0) { // if our loop did not terminate cause of EOF process as first read
						stack_insert(&stack, read);
						if (!is_unprocessable(read)) {
							// put read in second
							second = read;
						}
						read = bam_init1();
					} else { break; }
				}
				else {
					stack_insert(&stack, read);
					read = bam_init1();
				}
			}
		} else {
			// We don't have a pair
			if (!unprocessable) {
				// put read in second
				second = read;
			}
			stack_insert(&stack, read);
			read = bam_init1();
		}
	}
	// TODO: process as single read here
	flush_stack_to_disk(&stack, state->fout, state->hout, false);
	
	bam_destroy1(read);
	return success;
}