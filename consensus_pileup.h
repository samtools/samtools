/*  consensus_pileup.h -- Pileup orientated data per consensus column

    Copyright (C) 2013-2016, 2020-2021 Genome Research Ltd.

    Author: James Bonfied <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <htslib/sam.h>

typedef struct pileup {
    // commonly used things together, to fit in a cache line (64 bytes)
    struct pileup *next;  // A link list, for active seqs
    void *cd;             // General purpose per-seq client-data
    int  eof;             // True if this sequence has finished
    int  qual;            // Current qual (for active seq only)
    char start;           // True if this is a new sequence
    char base;            // Current base (for active seq only) in ASCII
    char ref_skip;        // True if the cause of eof or start is cigar N
    char padding;         // True if the base was added due to another seq
    int  base4;           // Base in 4-bit notation (0-15)
    hts_pos_t pos;        // Current unpadded position in seq
    int  nth;             // nth base at unpadded position 'pos'
    int  b_is_rev;        // 0 => fwd, 1 => rev
    int  seq_offset;      // Current base position in s->seq[] array.

    unsigned char *b_qual;// cached bam_qual
    unsigned char *b_seq; // cached bam_seq

    // --- 64 bytes
    struct pileup *eofn;  // p->eof set, next eof member
    struct pileup *eofl;  // last non-eof that points to p with p->eof

    uint32_t *b_cigar;    // cached bam_cigar

    int  cigar_ind;       // Current location in s->alignment cigar str
    int  cigar_op;        // Current cigar operation
    int  cigar_len;       // Remaining length of this cigar op

    int  first_del;       // Used when first base is a deletion

    bam1_t b;             // Bam entry associated with struct
} pileup_t;

/*
 * The pileup loop executes and calls callbacks to perform the work.
 *
 * seq_fetch returns the next sequence.  Return 0 from this indicates no
 *   more data.
 *
 * seq_init is called, if non-NULL, when a sequence is added to the pileup,
 * seq_free likewise, if non-NULL, is called when a sequence is removed
 *   from the pileup.
 * These two functions are akin to the constructor and destructors added
 * to mpileup.
 *
 * seq_column is the primary work horse which is executed for each
 *   reference position, and for each inserted base per ref pos.
 *
 * If we were to invert this from a loop generating callbacks to a polled
 * style interface like mpileup, then the seq_column bit would be dropped
 * and replaced by the returned pileup and associated parameters.
 */
int pileup_loop(samFile *fp,
                sam_hdr_t *h,
                int (*seq_fetch)(void *client_data,
                                 samFile *fp,
                                 sam_hdr_t *h,
                                 bam1_t *b),
                int (*seq_init)(void *client_data,
                                samFile *fp,
                                sam_hdr_t *h,
                                pileup_t *p),
                int (*seq_column)(void *client_data,
                                  samFile *fp,
                                  sam_hdr_t *h,
                                  pileup_t *p,
                                  int depth,
                                  hts_pos_t pos,
                                  int nth,
                                  int is_insert),
                void (*seq_free)(void *client_data,
                                 samFile *fp,
                                 sam_hdr_t *h,
                                 pileup_t *p),
                void *client_data);
