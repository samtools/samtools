/*  bam_markdup.c -- Mark duplicates from a coord sorted file that has gone
                     through fixmates with the mate scoring option on.

    Copyright (C) 2017 Genome Research Ltd.

    Author: Andrew Whitwham <aw7@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE
*/

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include <ctype.h>
#include "htslib/thread_pool.h"
#include "htslib/sam.h"
#include "sam_opts.h"
#include "samtools.h"
#include "htslib/khash.h"
#include "htslib/klist.h"

typedef struct {
    int32_t single;
    int32_t this_ref;
    int32_t this_coord;
    int32_t other_ref;
    int32_t other_coord;
    int32_t leftmost;
    int32_t orientation;
} key_data_t;

typedef struct {
    bam1_t *p;
} in_hash_t;

typedef struct {
    bam1_t *b;
    int32_t pos;
    key_data_t pair_key;
    key_data_t single_key;
} read_queue_t;



static khint32_t do_hash(unsigned char *key, khint32_t len);

static khint_t hash_key(key_data_t key) {
    int i = 0;
    khint_t hash;

    if (key.single) {
        unsigned char sig[12];

        memcpy(sig + i, &key.this_ref, 4);      i += 4;
        memcpy(sig + i, &key.this_coord, 4);    i += 4;
        memcpy(sig + i, &key.orientation, 4);   i += 4;

        hash = do_hash(sig, i);
    } else {
        unsigned char sig[24];

        memcpy(sig + i, &key.this_ref, 4);      i += 4;
        memcpy(sig + i, &key.this_coord, 4);    i += 4;
        memcpy(sig + i, &key.other_ref, 4);     i += 4;
        memcpy(sig + i, &key.other_coord, 4);   i += 4;
        memcpy(sig + i, &key.leftmost, 4);      i += 4;
        memcpy(sig + i, &key.orientation, 4);   i += 4;

        hash = do_hash(sig, i);
    }

    return hash;
}


static int key_equal(key_data_t a, key_data_t b) {
    int match = 1;

    if (a.this_coord != b.this_coord)
        match = 0;
    else if (a.orientation != b.orientation)
        match = 0;
    else if (a.this_ref != b.this_ref)
        match = 0;
    else if (a.single != b.single)
        match = 0;

    if (!a.single) {
        if (a.other_coord != b.other_coord)
            match = 0;
        else if (a.leftmost != b.leftmost)
            match = 0;
        else if (a.other_ref != b.other_ref)
            match = 0;
    }

    return match;
}


#define __free_queue_element(p)
#define O_FF 2
#define O_RR 3
#define O_FR 5
#define O_RF 7

KHASH_INIT(reads, key_data_t, in_hash_t, 1, hash_key, key_equal) // read map hash
KLIST_INIT(read_queue, read_queue_t, __free_queue_element) // the reads buffer


/* Calculate the mate's unclipped start based on position and cigar string from MC tag. */

static int32_t unclipped_other_start(int32_t op, char *cigar) {
    char *c = cigar;
    int32_t clipped = 0;

    while (*c && *c != '*') {
        long num = 0;

        if (isdigit((int)*c)) {
            num = strtol(c, &c, 10);
        } else {
            num = 1;
        }

        if (*c == 'S' || *c == 'H') { // clips
            clipped += num;
        } else {
            break;
        }

        c++;
    }

    return op - clipped + 1;
}


/* Calculate the current read's start based on the stored cigar string. */

static int32_t unclipped_start(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int32_t clipped = 0;
    uint32_t i;

    for (i = 0; i < b->core.n_cigar; i++) {
        char c = bam_cigar_opchr(cigar[i]);

        if (c == 'S' || c == 'H') { // clips
            clipped += bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }

    return b->core.pos - clipped + 1;
}


/* Calculate the mate's unclipped end based on start position and cigar string from MC tag.*/

static int32_t unclipped_other_end(int32_t op, char *cigar) {
    char *c = cigar;
    int32_t refpos = 0;
    int skip = 1;

    while (*c && *c != '*') {
        long num = 0;

        if (isdigit((int)*c)) {
            num = strtol(c, &c, 10);
        } else {
            num = 1;
        }

        switch (*c) {
            case 'M':
            case 'D':
            case 'N':
            case '=':
            case 'X':
                refpos += num;
                skip = 0; // ignore initial clips
            break;

            case 'S':
            case 'H':
                if (!skip) {
                refpos += num;
            }
            break;
        }

        c++;
   }

    return  op + refpos;
}


/* Calculate the current read's end based on the stored cigar string. */

static int32_t unclipped_end(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int32_t end_pos, clipped = 0;
    int32_t i;

    end_pos = bam_endpos(b);

    // now get the clipped end bases (if any)
    // if we get to the beginning of the cigar string
    // without hitting a non-clip then the results are meaningless
    for (i = b->core.n_cigar - 1; i >= 0; i--) {
        char c = bam_cigar_opchr(cigar[i]);

        if (c == 'S' || c == 'H') { // clips
            clipped += bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }

    return end_pos + clipped;
}


/* The Bob Jenkins one_at_a_time hash to reduce the key to a 32 bit value. */

static khint32_t do_hash(unsigned char *key, khint32_t len) {
    khint32_t   hash, i;

    for (hash = 0, i = 0; i < len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }

    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

    return hash;
}


/* Get mate score from tag. */

static int64_t get_mate_score(bam1_t *b) {
    uint8_t *data;
    int64_t score;

    if ((data = bam_aux_get(b, "ms"))) {
        score = bam_aux2i(data);
    } else {
        fprintf(stderr, "[markdup] error: no ms score tag.\n");
        return -1;
    }

    return score;
}


/* Calc current score from quality. */

static int64_t calc_score(bam1_t *b)
{
    int64_t score = 0;
    uint8_t  *qual = bam_get_qual(b);
    int i;

    for (i = 0; i < b->core.l_qseq; i++) {
        if (qual[i] >= 15) score += qual[i];
    }

    return score;
}


/* Create a signature hash of the current read and its pair.
   Uses the unclipped start (or end depending on orientation),
   the reference id, orientation and whether the current
   read is leftmost of the pair. */

static int make_pair_key(key_data_t *key, bam1_t *bam) {
    int32_t this_ref, this_coord, this_end;
    int32_t other_ref, other_coord, other_end;
    int32_t orientation, leftmost;
    uint8_t *data;
    char *cig;

    this_ref    = bam->core.tid + 1; // avoid a 0 being put into the hash
    other_ref   = bam->core.mtid + 1;

    this_coord = unclipped_start(bam);
    this_end   = unclipped_end(bam);

    if ((data = bam_aux_get(bam, "MC"))) {
        cig = bam_aux2Z(data);
        other_end   = unclipped_other_end(bam->core.mpos, cig);
        other_coord = unclipped_other_start(bam->core.mpos, cig);
    } else {
        fprintf(stderr, "[markdup] error: no MC tag.\n");
        return 1;
    }

    // work out orientations
    if (this_ref != other_ref) {
        leftmost = this_ref < other_ref;
    } else {
        if (bam_is_rev(bam) == bam_is_mrev(bam)) {
            if (!bam_is_rev(bam)) {
                leftmost = this_coord <= other_coord;
            } else {
                leftmost = this_end <= other_end;
            }
        } else {
            if (bam_is_rev(bam)) {
                leftmost = this_end <= other_coord;
            } else {
                leftmost = this_coord <= other_end;
            }
        }
    }

    // pair orientation
    if (leftmost) {
        if (bam_is_rev(bam) == bam_is_mrev(bam)) {
            other_coord = other_end;

            if (!bam_is_rev(bam)) {
                if (bam->core.flag & BAM_FREAD1) {
                    orientation = O_FF;
                } else {
                    orientation = O_RR;
                }
            } else {
                if (bam->core.flag & BAM_FREAD1) {
                    orientation = O_RR;
                } else {
                    orientation = O_FF;
                }
            }
        } else {
            if (!bam_is_rev(bam)) {
                orientation = O_FR;
                other_coord = other_end;
            } else {
                orientation = O_RF;
                this_coord = this_end;
            }
        }
    } else {
        if (bam_is_rev(bam) == bam_is_mrev(bam)) {
            this_coord = this_end;

            if (!bam_is_rev(bam)) {
                if (bam->core.flag & BAM_FREAD1) {
                    orientation = O_RR;
                } else {
                    orientation = O_FF;
                }
            } else {
                if (bam->core.flag & BAM_FREAD1) {
                    orientation = O_FF;
                } else {
                    orientation = O_RR;
                }
            }
        } else {
            if (!bam_is_rev(bam)) {
                orientation = O_RF;
                other_coord = other_end;
            } else {
                orientation = O_FR;
                this_coord = this_end;
            }
        }
    }

    if (!leftmost)
        leftmost = 13;
    else
        leftmost = 11;

    key->single        = 0;
    key->this_ref      = this_ref;
    key->this_coord    = this_coord;
    key->other_ref     = other_ref;
    key->other_coord   = other_coord;
    key->leftmost      = leftmost;
    key->orientation   = orientation;

    return 0;
}


/* Create a signature hash of single read (or read with an unmatched pair).
   Uses unclipped start (or end depending on orientation), reference id,
   and orientation. */

static void make_single_key(key_data_t *key, bam1_t *bam) {
    int32_t this_ref, this_coord;
    int32_t orientation;

    this_ref = bam->core.tid + 1; // avoid a 0 being put into the hash

    if (bam_is_rev(bam)) {
        this_coord = unclipped_end(bam);
        orientation = O_RR;
    } else {
        this_coord = unclipped_start(bam);
        orientation = O_FF;
    }

    key->single        = 1;
    key->this_ref      = this_ref;
    key->this_coord    = this_coord;
    key->orientation   = orientation;
}


/* Compare the reads near each other (coordinate sorted) and try to spot the duplicates.
   Generally the highest quality scoring is chosen as the original and all others the duplicates.
   The score is based on the sum of the quality values (<= 15) of the read and its mate (if any).
   While single reads are compared to only one read of a pair, the pair will chosen as the original.
   The comparison is done on position and orientation, see above for details. */

static int bam_mark_duplicates(samFile *in, samFile *out, int remove_dups, int32_t max_length, int do_stats) {
    bam_hdr_t *header;
    khiter_t k;
    khash_t(reads) *pair_hash        = kh_init(reads);
    khash_t(reads) *single_hash      = kh_init(reads);
    klist_t(read_queue) *read_buffer = kl_init(read_queue);
    kliter_t(read_queue) *rq;
    int32_t prev_tid, prev_coord;
    read_queue_t *in_read;
    int ret;
    int reading, writing, excluded, duplicate, single, pair, single_dup, examined;

    if ((header = sam_hdr_read(in)) == NULL) {
        fprintf(stderr, "[markdup] error reading header\n");
        return 1;
    }

    // accept unknown, unsorted or coordinate sort order, but error on queryname sorted.
    // only really works on coordinate sorted files.
    if ((header->l_text > 3) && (strncmp(header->text, "@HD", 3) == 0)) {
        char *p, *q;

       p = strstr(header->text, "\tSO:queryname");
       q = strchr(header->text, '\n');

       // looking for SO:queryname within @HD only
       // (e.g. must ignore in a @CO comment line later in header)
       if ((p != 0) && (p < q)) {
           fprintf(stderr, "[markdup] error: queryname sorted, must be sorted by coordinate.\n");
           return 1;
       }
    }

    if (sam_hdr_write(out, header) < 0) {
        fprintf(stderr, "[markdup] error writing header.\n");
        return 1;
    }

    // used for coordinate order checks
    prev_tid = prev_coord = 0;

    // get the buffer going
    in_read = kl_pushp(read_queue, read_buffer);


    if ((in_read->b = bam_init1()) == NULL) {
        fprintf(stderr, "[markdup] error: unable to allocate memory for alignment.\n");
        return 1;
    }

    reading = writing = excluded = single_dup = duplicate = examined = pair = single = 0;

    while ((ret = sam_read1(in, header, in_read->b)) >= 0) {

        // do some basic coordinate order checks
        if (in_read->b->core.tid >= 0) { // -1 for unmapped reads
            if (in_read->b->core.tid < prev_tid ||
               ((in_read->b->core.tid == prev_tid) && (in_read->b->core.pos < prev_coord))) {
                fprintf(stderr, "[markdup] error: bad coordinate order.\n");
                return 1;
            }
        }

        prev_coord = in_read->pos = in_read->b->core.pos;
        prev_tid   =  in_read->b->core.tid;
        in_read->pair_key.single   = 1;
        in_read->single_key.single = 0;

        reading++;

        // read must not be secondary, supplementary, unmapped or failed QC
        if (!(in_read->b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FQCFAIL))) {
            examined++;

            // look at the pairs first
            if ((in_read->b->core.flag & BAM_FPAIRED) && !(in_read->b->core.flag & BAM_FMUNMAP)) {
                int ret, mate_tmp;
                key_data_t pair_key;
                key_data_t single_key;
                in_hash_t *bp;

                if (make_pair_key(&pair_key, in_read->b)) {
                    fprintf(stderr, "[markdup] error: unable to assign pair hash key.\n");
                    return 1;
                }

                make_single_key(&single_key, in_read->b);

                pair++;
                in_read->pos = single_key.this_coord; // cigar/orientation modified pos

                // put in singles hash for checking against non paired reads
                k = kh_put(reads, single_hash, single_key, &ret);

                if (ret > 0) { // new
                    // add to single duplicate hash
                    bp = &kh_val(single_hash, k);
                    bp->p = in_read->b;
                    in_read->single_key = single_key;
                } else if (ret == 0) { // exists
                    // look at singles only for duplication marking
                    bp = &kh_val(single_hash, k);

                    if (!(bp->p->core.flag & BAM_FPAIRED) || (bp->p->core.flag & BAM_FMUNMAP)) {
                        bam1_t *dup = bp->p;

                        // singleton will always be marked duplicate even if
                        // scores more than one read of the pair

                        bp->p = in_read->b;
                        dup->core.flag |= BAM_FDUP;
                        single_dup++;
                    }
                } else {
                    fprintf(stderr, "[markdup] error: single hashing failure.\n");
                    return 1;
                }

                // now do the pair
                k = kh_put(reads, pair_hash, pair_key, &ret);

                if (ret > 0) { // new
                    // add to the pair hash
                    bp = &kh_val(pair_hash, k);
                    bp->p = in_read->b;
                    in_read->pair_key = pair_key;
                } else if (ret == 0) {
                    int64_t old_score, new_score, tie_add = 0;
                    bam1_t *dup;

                    bp = &kh_val(pair_hash, k);

                    if ((mate_tmp = get_mate_score(bp->p)) == -1) {
                        fprintf(stderr, "[markdup] error: no ms score tag.\n");
                        return 1;
                    } else {
                        old_score = calc_score(bp->p) + mate_tmp;
                    }

                    if ((mate_tmp = get_mate_score(in_read->b)) == -1) {
                        fprintf(stderr, "[markdup] error: no ms score tag.\n");
                        return 1;
                    } else {
                        new_score = calc_score(in_read->b) + mate_tmp;
                    }

                    // choose the highest score as the original
                    // and add it to the pair hash, mark the other as duplicate

                    if (new_score == old_score) {
                        if (strcmp(bam_get_qname(in_read->b), bam_get_qname(bp->p)) < 0) {
                            tie_add = 1;
                        } else {
                            tie_add = -1;
                        }
                    }

                    if (new_score + tie_add > old_score) { // swap reads
                        dup = bp->p;
                        bp->p = in_read->b;
                    } else {
                        dup = in_read->b;
                    }

                    dup->core.flag |= BAM_FDUP;

                    duplicate++;
                } else {
                    fprintf(stderr, "[markdup] error: pair hashing failure.\n");
                    return 1;
                }
            } else { // do the single (or effectively single) reads
                int ret;
                key_data_t single_key;
                in_hash_t *bp;

                make_single_key(&single_key, in_read->b);

                single++;
                in_read->pos = single_key.this_coord; // cigar/orientation modified pos

                k = kh_put(reads, single_hash, single_key, &ret);

                if (ret > 0) { // new
                    bp = &kh_val(single_hash, k);
                    bp->p = in_read->b;
                    in_read->single_key = single_key;
                } else if (ret == 0) { // exists
                    bp = &kh_val(single_hash, k);

                    if ((bp->p->core.flag & BAM_FPAIRED) && !(bp->p->core.flag & BAM_FMUNMAP)) {
                        // if matched against one of a pair just mark as duplicate
                        in_read->b->core.flag |= BAM_FDUP;
                    } else {
                        int64_t old_score, new_score;
                        bam1_t *dup;

                        old_score = calc_score(bp->p);
                        new_score = calc_score(in_read->b);

                        // choose the highest score as the original, add it
                        // to the single hash and mark the other as duplicate
                        if (new_score > old_score) { // swap reads
                            dup = bp->p;
                            bp->p = in_read->b;
                        } else {
                            dup = in_read->b;
                        }

                        dup->core.flag |= BAM_FDUP;
                    }

                    single_dup++;
                } else {
                    fprintf(stderr, "[markdup] error: single hashing failure.\n");
                    return 1;
                }
            }
        } else {
            excluded++;
        }

        // loop through the stored reads and write out those we
        // no longer need
        rq = kl_begin(read_buffer);
        while (rq != kl_end(read_buffer)) {
            in_read = &kl_val(rq);

            /* keep a moving window of reads based on coordinates and max read length.  Any unaligned reads
               should just be written as they cannot be matched as duplicates. */
            if (in_read->pos + max_length > prev_coord && in_read->b->core.tid == prev_tid && (prev_tid != -1 || prev_coord != -1)) {
                break;
            }

            if (!remove_dups || !(in_read->b->core.flag & BAM_FDUP)) {
                if (sam_write1(out, header, in_read->b) < 0) {
                    fprintf(stderr, "[markdup] error: writing output failed.\n");
                    return 1;
                }

                writing++;
            }

            // remove from hash
            if (in_read->pair_key.single == 0) {
                k = kh_get(reads, pair_hash, in_read->pair_key);
                kh_del(reads, pair_hash, k);
            }

            if (in_read->single_key.single == 1) {
                k = kh_get(reads, single_hash, in_read->single_key);
                kh_del(reads, single_hash, k);
            }

            kl_shift(read_queue, read_buffer, NULL);
            bam_destroy1(in_read->b);
            rq = kl_begin(read_buffer);
        }

        // set the next one up for reading
        in_read = kl_pushp(read_queue, read_buffer);

        if ((in_read->b = bam_init1()) == NULL) {
            fprintf(stderr, "[markdup] error: unable to allocate memory for alignment.\n");
            return 1;
        }
    }

    if (ret < -1) {
        fprintf(stderr, "[markdup] error: truncated input file.\n");
        return 1;
    }

    // write out the end of the list
    rq = kl_begin(read_buffer);
    while (rq != kl_end(read_buffer)) {
        in_read = &kl_val(rq);

        if (bam_get_qname(in_read->b)) { // last entry will be blank
            if (!remove_dups || !(in_read->b->core.flag & BAM_FDUP)) {
                if (sam_write1(out, header, in_read->b) < 0) {
                    fprintf(stderr, "[markdup] error: writing final output failed.\n");
                    return 1;
                }

                writing++;
            }
        }

        kl_shift(read_queue, read_buffer, NULL);
        bam_destroy1(in_read->b);
        rq = kl_begin(read_buffer);
    }

    if (do_stats) {
        fprintf(stderr, "READ %d WRITTEN %d \n"
            "EXCLUDED %d EXAMINED %d\n"
            "PAIRED %d SINGLE %d\n"
            "DULPICATE PAIR %d DUPLICATE SINGLE %d\n"
            "DUPLICATE TOTAL %d\n", reading, writing, excluded, examined, pair, single,
                                duplicate, single_dup, single_dup + duplicate);
    }

    kh_destroy(reads, pair_hash);
    kh_destroy(reads, single_hash);
    kl_destroy(read_queue, read_buffer);
    bam_hdr_destroy(header);

    return 0;
}


static int markdup_usage(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  samtools markdup <input.bam> <output.bam>\n\n");
    fprintf(stderr, "Option: \n");
    fprintf(stderr, "  -r           Remove duplicate reads\n");
    fprintf(stderr, "  -l           Max read length (default 300 bases)\n");
    fprintf(stderr, "  -s           Report stats.\n");

    sam_global_opt_help(stderr, "-.O..@");

    fprintf(stderr, "\nThe input file must be coordinate sorted and must have gone"
                     " through fixmates with the mate scoring option on.\n");

    return 1;
}


int bam_markdup(int argc, char **argv) {
    int c, ret, remove_dups = 0, report_stats = 0;
    int32_t max_length = 300;
    samFile *in = NULL, *out = NULL;
    char wmode[3] = {'w', 'b', 0};
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    htsThreadPool p = {NULL, 0};

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "rsl:O:@:", lopts, NULL)) >= 0) {
        switch (c) {
            case 'r': remove_dups = 1; break;
            case 'l': max_length = atoi(optarg); break;
            case 's': report_stats = 1; break;
            default: if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
            case '?': return markdup_usage();
        }
    }

    if (optind + 2 > argc)
        return markdup_usage();

    in = sam_open_format(argv[optind], "r", &ga.in);

    if (!in) {
        print_error_errno("markdup", "failed to open \"%s\" for input", argv[optind]);
        return 1;
    }

    sam_open_mode(wmode + 1, argv[optind + 1], NULL);
    out = sam_open_format(argv[optind + 1], wmode, &ga.out);

    if (!out) {
        print_error_errno("markdup", "failed to open \"%s\" for output", argv[optind + 1]);
        return 1;
    }

    if (ga.nthreads > 0)  {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "[markdup] error creating thread pool\n");
            return 1;
        }

        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    }

    // actual stuff happens here
    ret = bam_mark_duplicates(in, out, remove_dups, max_length, report_stats);

    sam_close(in);

    if (sam_close(out) < 0) {
        fprintf(stderr, "[markdup] error closing output file\n");
        ret = 1;
    }

    if (p.pool) hts_tpool_destroy(p.pool);

    sam_global_args_free(&ga);

    return ret;
}
