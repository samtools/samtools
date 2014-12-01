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
#include <htslib/klist.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdio.h>
#include "bam_mkdup.h"
#include "pos_buffer.h"


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

typedef struct namescore {
	char* name;
	int score;
} namescore_t;

// FIXME: write an alternative for platforms that don't have a 128 bit intrinsic
// FIXME: this is the 64bit hash find a good one for 128 bit numbers
#define possig_hash_func(key) (khint32_t)((key.bits)>>33^(key.bits)^(key.bits)<<11)

#define possig_hash_equal(a, b) ((a.bits) == (b.bits))

KHASH_INIT(signame, possig_t, namescore_t, 1, possig_hash_func, possig_hash_equal)

typedef struct endscore {
	possig_part_t possig_part;
	int score;
} endscore_t;

KHASH_MAP_INIT_STR(nameqs, endscore_t)
KHASH_MAP_INIT_STR(name, char)

#define __free_bam1_t(p)
KLIST_INIT(read, bam1_t*,__free_bam1_t)

static inline bool bam_to_possig_part(bam1_t* record, possig_part_t* part)
{
	if (record->core.tid > INT32_MAX) return false;
	part->orient = ((record->core.flag&BAM_FREVERSE) == BAM_FREVERSE);
	part->tid = record->core.tid;
	return true;
}

static inline bool part_bam_to_possig(possig_part_t* first_pos, bam1_t* second_pos, possig_t* sig)
{
	// Get 5 prime coord of each part of reads
	int32_t spos = (second_pos->core.flag&BAM_FREVERSE) ? bam_endpos(second_pos) : second_pos->core.pos;
	
	// Make the position signature with the first coordinate read in the first position
	if (first_pos->tid == second_pos->core.tid) {
		if (first_pos->pos < spos) {
			memcpy(&sig->field.first, first_pos, sizeof(possig_part_t));
			sig->field.second.pos = spos;
			return (bam_to_possig_part(second_pos, &sig->field.second));
		} else {
			memcpy(&sig->field.second, first_pos, sizeof(possig_part_t));
			sig->field.first.pos = spos;
			return (bam_to_possig_part(second_pos, &sig->field.first));
		}
	} else {
		if (first_pos->tid < second_pos->core.tid) {
			memcpy(&sig->field.first, first_pos, sizeof(possig_part_t));
			sig->field.second.pos = spos;
			return (bam_to_possig_part(second_pos, &sig->field.second));
		} else {
			memcpy(&sig->field.second, first_pos, sizeof(possig_part_t));
			sig->field.first.pos = spos;
			return (bam_to_possig_part(second_pos, &sig->field.first));
		}
	}
}

/*
 * To process reads in one location
 *
 * TODO
 *
 * 2nd pass
 * for each read r in b
 * 	if r.qname in kill_hash
 * 		mark r
 * 	endif
 * 	write r
 * end for
 */

// Returns total of all quality values in read
static inline int sum_qual(const bam1_t *b)
{
	int i, q;
	uint8_t *qual = bam_get_qual(b);
	for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
	return q;
}


bool process_coordsorted(/* HACK:const*/ state_t* state, const char* BIG_DIRTY_HACK)
{
	khash_t(signame)* ps_hash = kh_init(signame);
	khash_t(nameqs)* name_hash = kh_init(nameqs);
	khash_t(name)* kill_hash = kh_init(name);
	bam1_t* read_first = bam_init1();
	bool success = true;
	size_t distinct_pos = 0;
	size_t marked = 0;
	size_t killed = 0;

	// First pass
	while (sam_read1(state->fin, state->hin, read_first) >= 0) {
		if (is_unprocessable(read_first)) continue; // Skip non-primaries
		
		// Have we seen the other half of this read?
		khiter_t k = kh_get(nameqs, name_hash, bam_get_qname(read_first));
		if ( k != kh_end(name_hash) ) {
			// It's already there
			int score = kh_value(name_hash, k).score + sum_qual(read_first);
			possig_t possig;
			if (!part_bam_to_possig(&kh_value(name_hash, k).possig_part, read_first, &possig)) {return false;}
			kh_del(nameqs, name_hash, k);

			// Have we seen this position pair before?
			khiter_t kp = kh_get(signame, ps_hash, possig );
			if ( kp != kh_end(ps_hash) ) {
				// Already there so get it
				char* to_kill;
				int kill_score;
				++marked;
				// Does what is in there have a better score than us?
				if (kh_value(ps_hash, kp).score < score) {
					// Kill what's already there and replace it
					to_kill = kh_value(ps_hash, kp).name;
					kill_score = kh_value(ps_hash, kp).score;

					kh_value(ps_hash, kp).name = strdup(bam_get_qname(read_first));
					kh_value(ps_hash, kp).score = score;
				} else {
					to_kill = strdup(bam_get_qname(read_first));
					kill_score = score;
					++killed;
				}
				// If it's not already in the kill hash kill it
				khiter_t kill_iter = kh_get(name, kill_hash, to_kill);
				if (kill_iter == kh_end(kill_hash)) {
					int dummy;
					kill_iter = kh_put(name, kill_hash, to_kill, &dummy );
				} else {free(to_kill);}
				kh_value(kill_hash, kill_iter) = kill_score;
			} else {
				int not_found_p;
				khiter_t kp = kh_put(signame, ps_hash, possig, &not_found_p );
				kh_value(ps_hash, kp).name = strdup(bam_get_qname(read_first));
				kh_value(ps_hash, kp).score = score;

				++distinct_pos;
			}
		} else {
			// Not there so make it
			int not_found;
			k = kh_put(nameqs, name_hash, strdup(bam_get_qname(read_first)), &not_found );
			// it's been inserted
			bam_to_possig_part(read_first, &kh_value(name_hash, k).possig_part);
			kh_value(name_hash, k).possig_part.pos = (read_first->core.flag&BAM_FREVERSE) ? bam_endpos(read_first) : read_first->core.pos;
			kh_value(name_hash, k).score = sum_qual(read_first);
		}
	}
	bam_destroy1(read_first);
	
	fprintf(stderr, "distinct_pos: %zu marked: %zu killed: %zu Number to kill size:%d buckets:%d\n", distinct_pos, marked, killed, kh_size(kill_hash), kh_n_buckets(kill_hash));
	
	// 2nd pass
	// HACK: PUT REWIND API HERE
	sam_close(state->fin);
	state->fin = sam_open(BIG_DIRTY_HACK, "rb");
	state->hin = sam_hdr_read(state->fin);
	// END/HACK
	bam1_t* read = bam_init1();
	while (sam_read1(state->fin, state->hin, read) >= 0) {
		if ( kh_end(kill_hash) != kh_get(name, kill_hash, bam_get_qname(read))) {
			read->core.flag |= BAM_FDUP; // Set the duplicate flag
		} else {
			read->core.flag &= ~BAM_FDUP;  // Strip any existing duplicate flags
		}
		sam_write1(state->fout, state->hout, read);
	}
	bam_destroy1(read);
	kh_destroy(name, kill_hash);
	kh_destroy(nameqs, name_hash);

	return success;
}
