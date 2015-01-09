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

typedef struct endscore {
	read_vector_t read_vector;
	int score;
} endscore_t;

KHASH_MAP_INIT_STR(nameqs, endscore_t)
KHASH_MAP_INIT_STR(name, char)

typedef struct fragname {
    char* name;
    int score;
} fragname_t;

#define read_vector_hash_func(key) (khint32_t)((*(uint64_t*)(&key))>>33^(*(uint64_t*)(&key))^(*(uint64_t*)(&key))<<11)
#define read_vector_hash_equal(a, b) (*(uint64_t*)(&a) == *(uint64_t*)(&b))


KHASH_INIT(fragscore, read_vector_t, fragname_t, 1, read_vector_hash_func, read_vector_hash_equal)

#define __free_bam1_t(p)
KLIST_INIT(read, bam1_t*,__free_bam1_t)


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
	pos_buffer_t* buf = pos_buffer_init();
	khash_t(nameqs)* name_hash = kh_init(nameqs); // Hash of score and positions indexed by name for pairing
	khash_t(name)* kill_hash = kh_init(name); // Hash of reads to kill by name
    kh_fragscore_t* bestfrag = kh_init(fragscore); // Best scoring fragment at position
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
			// We have so now we have both parts of read
			int score = kh_value(name_hash, k).score + sum_qual(read_first);
			read_vector_t curr;
			if (!bam_to_read_vector(read_first, &curr)) abort();

            // put both halves into frag with INT_MAX score to mark non-paired fragments at that position
            int ret;
            khiter_t val;
            val = kh_put(fragscore, bestfrag, kh_value(name_hash, k).read_vector, &ret);
            kh_value(bestfrag, val).name = NULL; // Don't need to know name
            kh_value(bestfrag, val).score = INT_MAX; // INT_MAX score
            val = kh_put(fragscore, bestfrag, curr, &ret);
            kh_value(bestfrag, val).name = NULL; // Don't need to know name
            kh_value(bestfrag, val).score = INT_MAX; // INT_MAX score

			// Have we seen this position pair before?
			char* kill = pos_buffer_insert(buf, kh_value(name_hash, k).read_vector, curr, score, bam_get_qname(read_first), read_first->core.tid, read_first->core.pos);
			kh_del(nameqs, name_hash, k);

			if ( kill != NULL ) {
				// Already there so get it
				int kill_score = 1;
				++marked;
				// If it's not already in the kill hash kill it
				khiter_t kill_iter = kh_get(name, kill_hash, kill);
				if (kill_iter == kh_end(kill_hash)) {
					int dummy;
					kill_iter = kh_put(name, kill_hash, kill, &dummy );
				} else {free(kill);}
				kh_value(kill_hash, kill_iter) = kill_score;
			} else {
				++distinct_pos;
			}
		} else {
			// Not there so make it
			int not_found;
			k = kh_put(nameqs, name_hash, strdup(bam_get_qname(read_first)), &not_found );
			// it's been inserted
			bam_to_read_vector(read_first, &kh_value(name_hash, k).read_vector);
			kh_value(name_hash, k).read_vector.pos = (read_first->core.flag&BAM_FREVERSE) ? bam_endpos(read_first) : read_first->core.pos;
			kh_value(name_hash, k).score = sum_qual(read_first);
		}
	}
	bam_destroy1(read_first);
    
    // now handle what's left in nameqs which is fragments
    khiter_t k = kh_begin(name_hash);
    while (k != kh_end(name_hash)) {
        if (!kh_exist(name_hash, k)) { ++k; continue; }
        khiter_t val = kh_get(fragscore, bestfrag, kh_value(name_hash, k).read_vector);
        if (val != kh_end(bestfrag)) {
            // Already in there somebody is going to die
            // Do we score better than it?
            char *kill;
            if ( kh_value(name_hash, k).score > kh_value(bestfrag, val).score) {
                kill = strdup(kh_value(bestfrag, val).name);
                // Update best score in bestfrag
                kh_value(bestfrag, val).name = strdup(kh_key(name_hash, k));
                kh_value(bestfrag, val).score = kh_value(name_hash, k).score;
            } else {
                kill = strdup(kh_key(name_hash, k));
            }
            // kill the worst scorer
            int kill_score = 1;
            ++killed;
            // If it's not already in the kill hash kill it
            khiter_t kill_iter = kh_get(name, kill_hash, kill);
            if (kill_iter == kh_end(kill_hash)) {
                int dummy;
                kill_iter = kh_put(name, kill_hash, kill, &dummy );
            } else {free(kill);}
            kh_value(kill_hash, kill_iter) = kill_score;
        } else {
            int ret = 0;
            val = kh_put(fragscore, bestfrag, kh_value(name_hash, k).read_vector, &ret);
            // Not seen this location before so add it to the frag hash
            kh_value(bestfrag, val).name = strdup(kh_key(name_hash, k));
            kh_value(bestfrag, val).score = kh_value(name_hash, k).score;
        }

        ++k;
    }
    // Deallocate memory in bestfrag and destroy it
    fragname_t fg;
    kh_foreach_value(bestfrag, fg, free(fg.name));
    kh_destroy(fragscore, bestfrag);
	
	fprintf(stderr, "distinct_pos: %zu marked pairs: %zu marked frags: %zu Number to kill size:%d buckets:%d\n", distinct_pos, marked, killed, kh_size(kill_hash), kh_n_buckets(kill_hash));
	
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
