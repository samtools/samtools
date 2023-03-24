/*  consensus__pileup.h -- Pileup orientated data per consensus column

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

#include <config.h>
#include <htslib/sam.h>

#ifdef __SSE__
#   include <xmmintrin.h>
#else
#   define _mm_prefetch(a,b)
#endif

#include "consensus_pileup.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define bam_strand(b)  (((b)->core.flag & BAM_FREVERSE) != 0)

/*
 * START_WITH_DEL is the mode that Gap5 uses when building this. It prepends
 * all cigar strings with 1D and decrements the position by one. (And then
 * has code to reverse this operation in the pileup handler.)
 *
 * The reason for this is that it means reads starting with an insertion work.
 * Otherwise the inserted bases are silently lost. (Try it with "samtools
 * mpileup" and you can see it has the same issue.)
 *
 * However it's probably not want most people expect.
 */
//#define START_WITH_DEL

/* --------------------------------------------------------------------------
 * The pileup code itself.
 *
 * This consists of the external pileup_loop() function, which takes a
 * sam/bam samfile_t pointer and a callback function. The callback function
 * is called once per column of aligned data (so once per base in an
 * insertion).
 *
 * Current known issues.
 * 1) zero length matches, ie 2S2S cause failures.
 * 2) Insertions at starts of sequences get included in the soft clip, so
 *    2S2I2M is treated as if it's 4S2M
 * 3) From 1 and 2 above, 1S1I2S becomes 2S2S which fails.
 */


/*
 * Fetches the next base => the nth base at unpadded position pos. (Nth can
 * be greater than 0 if we have an insertion in this column). Do not call this
 * with pos/nth lower than the previous query, although higher is better.
 * (This allows it to be initialised at base 0.)
 *
 * Stores the result in base and also updates is_insert to indicate that
 * this sequence still has more bases in this position beyond the current
 * nth parameter.
 *
 * Returns 1 if a base was fetched
 *         0 if not (eg ran off the end of sequence)
 */
static int get_next_base(pileup_t *p, hts_pos_t pos, int nth, int *is_insert) {
    bam1_t *b = &p->b;
    int op = p->cigar_op;

    p->start -= p->start>0;
    if (p->first_del && op != BAM_CPAD)
        p->first_del = 0;

    *is_insert = 0;

    /* Find pos first */
    while (p->pos < pos) {
        p->nth = 0;

        if (p->cigar_len == 0) {
            if (p->cigar_ind >= b->core.n_cigar) {
                p->eof = 1;
                return 0;
            }

            op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
            p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
            p->cigar_ind++;
        }

        if ((op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
            && p->cigar_len <= pos - p->pos) {
            p->seq_offset += p->cigar_len;
            p->pos += p->cigar_len;
            p->cigar_len = 0;
        } else {
            switch (op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                p->seq_offset++;
                /* Fall through */
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                p->pos++;
                p->cigar_len--;
                break;

            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                p->seq_offset += p->cigar_len;
                /* Fall through */
            case BAM_CPAD:
            case BAM_CHARD_CLIP:
                p->cigar_len = 0;
                break;

            default:
                fprintf(stderr, "Unhandled cigar_op %d\n", op);
                return -1;
            }
        }
    }

    /* Now at pos, find nth base */
    while (p->nth < nth) {
        if (p->cigar_len == 0) {
            if (p->cigar_ind >= b->core.n_cigar) {
                p->eof = 1;
                return 0; /* off end of seq */
            }

            op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
            p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
            p->cigar_ind++;
        }

        switch (op) {
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
        case BAM_CSOFT_CLIP:
        case BAM_CDEL:
        case BAM_CREF_SKIP:
            goto at_nth; /* sorry, but it's fast! */

        case BAM_CINS:
            p->seq_offset++;
            /* Fall through */
        case BAM_CPAD:
            p->cigar_len--;
            p->nth++;
            break;

        case BAM_CHARD_CLIP:
            p->cigar_len = 0;
            break;

        default:
            fprintf(stderr, "Unhandled cigar_op %d\n", op);
            return -1;
        }
    }
 at_nth:

    /* Fill out base & qual fields */
    p->ref_skip = 0;
    if (p->nth < nth && op != BAM_CINS) {
        //p->base = '-';
        p->base = '*';
        p->base4 = 16;
        p->padding = 1;
        if (p->seq_offset < b->core.l_qseq)
            p->qual = MIN(p->qual, p->b_qual[p->seq_offset+1]);
        else
            p->qual = 0;
    } else {
        p->padding = 0;
        switch(op) {
        case BAM_CDEL:
            p->base = '*';
            p->base4 = 16;
            if (p->seq_offset+1 < b->core.l_qseq)
                p->qual = MIN(p->qual, p->b_qual[p->seq_offset+1]);
            else
                p->qual = MIN(p->qual, p->b_qual[p->seq_offset]);
            break;

        case BAM_CPAD:
            //p->base = '+';
            p->base = '*';
            p->base4 = 16;
            if (p->seq_offset+1 < b->core.l_qseq)
                p->qual = MIN(p->qual, p->b_qual[p->seq_offset+1]);
            else
                p->qual = MIN(p->qual, p->b_qual[p->seq_offset]);
            break;

        case BAM_CREF_SKIP:
            p->base = '.';
            p->base4 = 0;
            p->qual = 0;
            /* end of fragment, but not sequence */
            p->eof = p->eof ? 2 : 3;
            p->ref_skip = 1;
            break;

        default:
            if (p->seq_offset < b->core.l_qseq) {
                p->qual = p->b_qual[p->seq_offset];
                p->base4 = p->b_seq[p->seq_offset/2] >>
                    ((~p->seq_offset&1)<<2) & 0xf;
                p->base = "NACMGRSVTWYHKDBN"[p->base4];
            } else {
                p->base = 'N';
                p->base4 = 15;
                p->qual = 0xff;
            }

            break;
        }
    }

    /* Handle moving out of N (skip) into sequence again */
    if (p->eof && p->base != '.') {
        p->start = 1;
        p->ref_skip = 1;
        p->eof = 0;
    }

    /* Starting with an indel needs a minor fudge */
    if (p->start && p->cigar_op == BAM_CDEL) {
        p->first_del = 1;
    }

    /* Check if next op is an insertion of some sort */
    if (p->cigar_len == 0) {
        if (p->cigar_ind < b->core.n_cigar) {
            op=p->cigar_op  = p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK;
            p->cigar_len = p->b_cigar[p->cigar_ind] >> BAM_CIGAR_SHIFT;
            p->cigar_ind++;
            if (op == BAM_CREF_SKIP) {
                p->eof = 3;
                p->ref_skip = 1;
            }
        } else {
            p->eof = 1;
        }
    }

    switch (op) {
    case BAM_CPAD:
    case BAM_CINS:
        *is_insert = p->cigar_len;
        break;

    case BAM_CSOFT_CLIP:
        /* Last op 'S' => eof */
        p->eof = (p->cigar_ind == b->core.n_cigar ||
                  (p->cigar_ind+1 == b->core.n_cigar &&
                   (p->b_cigar[p->cigar_ind] & BAM_CIGAR_MASK)
                   == BAM_CHARD_CLIP))
            ? 1
            : 0;
        break;

    case BAM_CHARD_CLIP:
        p->eof = 1;
        break;

    default:
        break;
    }

    return 1;
}

/*
 * Loops through a set of supplied ranges producing columns of data.
 * When found, it calls func with clientdata as a callback. Func should
 * return 0 for success and non-zero for failure. seq_init() is called
 * on each new entry before we start processing it. It should return 0 or 1
 * to indicate reject or accept status (eg to filter unmapped data).
 * If seq_init() returns -1 we abort the pileup_loop with an error.
 * seq_init may be NULL.
 *
 * Returns 0 on success
 *        -1 on failure
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
                void *client_data) {
    int ret = -1;
    pileup_t *phead = NULL, *p, *pfree = NULL, *last, *next, *ptail = NULL;
    pileup_t *pnew = NULL;
    int is_insert, nth = 0, r;
    hts_pos_t col = 0;
    int last_ref = -1;

    /* FIXME: allow for start/stop boundaries rather than consuming all data */

    if (NULL == (pnew = calloc(1, sizeof(*p))))
        return -1;

    do {
        bam1_t *b;
        hts_pos_t pos;

        r = seq_fetch(client_data, fp, h, &pnew->b);
        if (r < -1) {
            fprintf(stderr, "bam_next_seq() failure.\n");
            goto error;
        }

        b = &pnew->b;

        /* Force realloc */
        //fp->bs = NULL;
        //fp->bs_size = 0;

        //r = samread(fp, pnew->b);
        if (r >= 0) {
            if (b->core.flag & BAM_FUNMAP)
                continue;

            if (b->core.tid == -1) {
                /* Another indicator for unmapped */
                continue;
            } else if (b->core.tid == last_ref) {
                pos = b->core.pos+1;
                //printf("New seq at pos %d @ %d %s\n", pos, b->core.tid,
                //       bam_name(b));
            } else {
                //printf("New ctg at pos %ld @ %d\n",b->core.pos+1,b->core.tid);
                pos = HTS_POS_MAX;
            }
        } else {
            pos = HTS_POS_MAX;
        }

        if (col > pos) {
            fprintf(stderr, "BAM/SAM file is not sorted by position. "
                    "Aborting\n");
            goto error;
        }

        /* Process data between the last column and our latest addition */
        while (col < pos && phead) {
            struct pileup *eof_head = NULL, *eofp = NULL;
            int v, ins, depth = 0;
            //printf("Col=%ld pos=%ld nth=%d\n", col, pos, nth);

            /* Pileup */
            is_insert = 0;
            pileup_t *pnext = phead ? phead->next : NULL;
            for (p = phead, last = NULL; p; p = pnext) {
#if 0
                // Simple prefetching
                pnext = p->next;
                if (pnext)
                    _mm_prefetch(pnext, _MM_HINT_T0);
#else
                // More complex prefetching => more instructions, but
                // usually faster.
                pnext = p->next;
                if (pnext) {
                    // start memory fetches; a big help on very deep data
                    if (pnext->next)
                        // struct 2 ahead
                        _mm_prefetch(pnext->next, _MM_HINT_T0);
                    // seq/qual 1 ahead
                    _mm_prefetch(pnext->b_qual + pnext->seq_offset,
                                 _MM_HINT_T0);
                    _mm_prefetch(pnext->b_seq  + pnext->seq_offset/2,
                                 _MM_HINT_T0);
                }
#endif

                if (!get_next_base(p, col, nth, &ins))
                    p->eof = 1;
                if (p->eof == 1) {
                    if (eofp)
                        eofp->eofn = p;
                    eofp = p;
                    eofp->eofl = last;
                    if (!eof_head)
                        eof_head = eofp;
                } else {
                    last = p;
                }

                if (is_insert < ins)
                    is_insert = ins;

                depth++;
            }
            if ((ptail = last) == NULL)
                ptail = phead;

            /* Call our function on phead linked list */
            v = seq_column(client_data, fp, h, phead, depth,
#ifdef START_WITH_DEL
                           col-1,
#else
                           col,
#endif
                           nth, is_insert);

            /* Remove dead seqs */
            for (p = eof_head ; p; p = p->eofn) {
                if (p->eofl)
                    p->eofl->next = p->next;
                else
                    phead = p->next;

                p->next = pfree;
                pfree = p;

                if (seq_free)
                    seq_free(client_data, fp, h, p);
            }

            if (v == 1)
                break; /* early abort */

            if (v != 0)
                goto error;

            /* Next column */
            if (is_insert) {
                nth++;
            } else {
                nth = 0;
                col++;
            }
        }

        /* May happen if we have a hole in the contig */
        col = pos;

        /* New contig */
        if (b && b->core.tid != last_ref) {
            last_ref = b->core.tid;
            pos = b->core.pos+1;
            nth = 0;
            col = pos;
        }

        /*
         * Add this seq.
         * Note: cigars starting with I or P ops (eg 2P3I10M) mean we have
         * alignment instructions that take place before the designated
         * starting location listed in the SAM file. They won't get included
         * in the callback function until they officially start, which is
         * already too late.
         *
         * So to workaround this, we prefix all CIGAR with 1D, move the
         * position by 1bp, and then force the callback code to remove
         * leaving pads (either P or D generated).
         *
         * Ie it's a level 10 hack!
         */
        if (r >= 0) {
            p = pnew;
            p->next       = NULL;
            p->cd         = NULL;
            p->eofn       = NULL;
            p->eofl       = NULL;
            p->start      = 2;
            p->eof        = 0;
#ifdef START_WITH_DEL
            p->pos        = pos-1;
            p->cigar_ind  = 0;
            p->b_cigar    = bam_get_cigar(&p->b);
            if ((p->b_cigar[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) {
                p->cigar_len  = p->b_cigar[0] >> BAM_CIGAR_SHIFT;
                p->cigar_op   = BAM_CHARD_CLIP;
                if ((p->b_cigar[1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
                    /* xHxS... => xHxS1D... */
                    p->b_cigar[0] = p->b_cigar[1];
                    p->b_cigar[1] = (1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
                } else {
                    /* xH... => xH1D... */
                    p->b_cigar[0] = (1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
                }
            } else {
                if ((p->b_cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
                    /* xS... => xS1D... */
                    p->cigar_len  = p->b_cigar[0] >> BAM_CIGAR_SHIFT;
                    p->cigar_op   = BAM_CSOFT_CLIP;
                    p->b_cigar[0] = (1 << BAM_CIGAR_SHIFT) | BAM_CDEL;
                } else {
                    /* ... => 1D... */
                    p->cigar_len  = 1;        /* was  0  */
                    p->cigar_op   = BAM_CDEL; /* was 'X' */
                }
            }
            p->seq_offset = -1;
            p->first_del  = 1;
#else
            p->pos        = pos-1;
            p->cigar_ind  = 0;
            p->b_cigar    = bam_get_cigar(&p->b);
            p->cigar_len  = 0;
            p->cigar_op   = -1;
            p->seq_offset = -1;
            p->first_del  = 0;
#endif
            p->b_is_rev   = bam_is_rev(&p->b);
            p->b_qual     = (uint8_t *)bam_get_qual(&p->b);
            p->b_seq      = (uint8_t *)bam_get_seq(&p->b);

            if (seq_init) {
                int v;
                v = seq_init(client_data, fp, h, p);
                if (v == -1)
                    goto error;

                if (v == 1) {
                    /* Keep this seq */
                    if (phead) {
                        ptail->next = p;
                    } else {
                        phead = p;
                    }
                    ptail = p;
                } else {
                    /* Push back on free list */
                    p->next = pfree;
                    pfree = p;
                }
            } else {
                if (phead)
                    ptail->next = p;
                else
                    phead = p;
                ptail = p;
            }

            /* Allocate the next pileup rec */
            if (pfree) {
                pnew = pfree;
                pfree = pfree->next;
            } else {
                if (NULL == (pnew = calloc(1, sizeof(*pnew))))
                    goto error;
            }
        }
    } while (r >= 0);

    ret = 0;
 error:

    if (pnew) {
        free(pnew->b.data);
        free(pnew);
    }

    /* Tidy up */
    for (p = pfree; p; p = next) {
        next = p->next;
        if (seq_free)
            seq_free(client_data, fp, h, p);
        free(p->b.data);
        free(p);
    }

    return ret;
}
