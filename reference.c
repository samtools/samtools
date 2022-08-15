/*  bam_reference.c -- extracts an embedded reference from a CRAM file,
                       or creates it from alignments plus MD:Z tags.

    Copyright (C) 2022 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include "htslib/sam.h"
#include "htslib/cram.h"
#include "samtools.h"
#include "sam_opts.h"


/*
 * There are two main modes of operation.
 *
 * 1. Extracting the reference from the CRAM file embed_ref blocks.
 * 2. Generation of reference by analysing consensus plus patches applied
 *    via MD tags.
 *
 * The first is very rapid, but only applies to a CRAM files generated with
 * the specific options (not commonly used and not the default).  The second
 * is a slow operation, but applies to any data type.
 *
 * This is also a testing ground for a future CRAM auto-embed-ref option that
 * permits the use of an embedded reference without having to first extract
 * the reference.  (Note this may require the creation of MD tags during
 * decode by use of an existing embedded reference, if the records don't
 * have an MD tag themselves, but that's an issue for htslib when we get
 * there.)
 */

/*
 * ---------------------------------------------------------------------------
 * Shared utility functions by both methods.
 */

#define haszero(x) (((x)-0x0101010101010101UL)&~(x)&0x8080808080808080UL)
#define MIN(a,b) ((a)<(b)?(a):(b))
static int dump_ref(sam_hdr_t *h, hts_itr_t *iter, int ref_id,
                    char *ref, uint64_t ref_len, FILE *fp, int verbose) {
    int N = 0;
    if (iter && iter->end >= HTS_POS_MAX)
        iter->end = ref_len;
    if (iter && (iter->beg > 0 || iter->end < ref_len)) {
        fprintf(fp, ">%s:%"PRIhts_pos"-%"PRIhts_pos"\n",
                sam_hdr_tid2name(h, ref_id), iter->beg+1, iter->end);
        ref += iter->beg;
        ref_len = MIN(ref_len, iter->end) - iter->beg;
    } else {
        fprintf(fp, ">%s\n", sam_hdr_tid2name(h, ref_id));
    }

    int i, j;
    uint64_t rem = ref_len;

    // Count coverage, purely for information purposes.
    // About 90% of dump_ref CPU is here, so maybe this isn't useful,
    // but this is still 3-4x faster than the obvious naive loop.
    //
    // Overall though it's only about 5% overhead of the entire process
    // (was ~20%).
    if (verbose) {
        int n4[8] = {0};
        for (j = 0; j < ref_len && (((uintptr_t) &ref[j] & 7) != 0); j++)
            N += ref[j] == 'N';
        uint64_t fast_end = ((ref_len - j) & ~7) + j;
        for (; j < fast_end; j+=8) {
            uint64_t i64 = *(uint64_t *)&ref[j];
            if (!haszero(i64 ^ 0x4e4e4e4e4e4e4e4eUL)) // 'N' <-> 0
                continue;

            n4[0] += ref[j+0] == 'N';
            n4[1] += ref[j+1] == 'N';
            n4[2] += ref[j+2] == 'N';
            n4[3] += ref[j+3] == 'N';
            n4[4] += ref[j+4] == 'N';
            n4[5] += ref[j+5] == 'N';
            n4[6] += ref[j+6] == 'N';
            n4[7] += ref[j+7] == 'N';
        }
        for (; j < ref_len; j++)
            N += ref[j] == 'N';
        N += n4[0]+n4[1]+n4[2]+n4[3]+
            n4[4]+n4[5]+n4[6]+n4[7];
    }

    // Format reference
    for (i = 0; i < ref_len; i += 60, rem -= 60) {
        int len = (int)(rem < 60 ? rem : 60);
        if (fwrite(ref, 1, len, fp) != len)
            return -1;
        putc('\n', fp);
        ref += 60;
    }

    if (verbose)
        fprintf(stderr, "Dump ref %d len %"PRId64", coverage %.2f%%\n",
                ref_id, ref_len, 100 - N*100.0 / ref_len);

    return 0;
}

/*
 * ---------------------------------------------------------------------------
 * CRAM embedded reference method of reference construction
 */

/*
 * Extracts an embedded reference from a sorted CRAM file.
 * Modelled on the CRAM container copy loop from bam_cat.c.
 */
static int cram2ref(samFile *in, sam_hdr_t *h, hts_idx_t *idx, char *reg,
                    FILE *outfp, int verbose) {
    cram_fd *in_c;
    cram_container *c = NULL;
    cram_block *blk = NULL;
    cram_block_slice_hdr *shdr = NULL;

    int curr_ref_id = -99;
    char *ref = NULL;
    uint64_t ref_len = 0;

    // We have no direct public API for seeking in CRAM to a specific
    // location by genome coordinates.  The sam_itr_query API is
    // designed for fetching records, rather than seeks to specific
    // file locations.
    //
    // TODO: consider exposing cram_range and cram_seek_to_refpos API.
    // After a sam_index_load which will add the index to infp, these
    // functions should seek direct to the start of a container.
    // Or use cram_index *e =cram_index_query(cram, tid, beg, NULL);
    //
    // However, fortuitously(?) sam_itr_querys calls cram_seek_to_refpos
    // so we can do a region query and let that do the initial seek.
    // We still need to do our own end-range detection though.

    hts_itr_t *iter = NULL;
    if (reg) {
        iter = sam_itr_querys(idx, h, reg);
        if (!iter) {
            print_error("reference", "failed to parse region '%s'", reg);
            goto err;
        }
    }

    in_c = in->fp.cram; // low level htslib abuse?
    int eor = 0;
    while (!eor && (c = cram_read_container(in_c))) {
        if (cram_container_is_empty(in_c)) {
            cram_block *blk;
            // Container compression header
            if (!(blk = cram_read_block(in_c)))
                goto err;
            cram_free_block(blk);
            cram_free_container(c);
            c = NULL; blk = NULL;
            continue;
        }

        // Container compression header; read and discard
        int32_t num_slices;
        if (!(blk = cram_read_block(in_c)))
            goto err;
        cram_free_block(blk);
        blk = NULL;

        // Container num_blocks can be invalid, due to a bug.
        // Instead we iterate in slice context instead.
        (void)cram_container_get_landmarks(c, &num_slices);
        int i, j;
        for (i = 0; i < num_slices; i++) {
            // Slice header
            if (!(blk = cram_read_block(in_c)))
                goto err;
            if (!(shdr = cram_decode_slice_header(in_c, blk)))
                goto err;
            cram_free_block(blk);
            blk = NULL;

            int num_blocks = cram_slice_hdr_get_num_blocks(shdr);
            int embed_id = cram_slice_hdr_get_embed_ref_id(shdr);
            int ref_id;
            hts_pos_t ref_start, ref_span;
            cram_slice_hdr_get_coords(shdr, &ref_id, &ref_start, &ref_span);

            if (iter) {
                if (iter->tid != ref_id || ref_start > iter->end) {
                    // Beyond end of specified region.
                    cram_free_slice_header(shdr);
                    eor = 1;
                    break;
                }
            }

            if (embed_id < 0 && ref_id != -1) {
                fprintf(stderr, "CRAM file has slice without embedded "
                        "reference\n");
                goto err;
            }

            if (ref_id != curr_ref_id) {
                if (curr_ref_id >= 0) {
                    if (dump_ref(h, iter, curr_ref_id, ref, ref_len,
                                 outfp, verbose) < 0)
                        goto err;
                }

                ref_len = sam_hdr_tid2len(h, ref_id);
                if (ref_len) {
                    char *ref2 = realloc(ref, ref_len);
                    if (!ref2)
                        goto err;
                    else
                        ref = ref2;
                    memset(ref, 'N', ref_len);
                }
                curr_ref_id = ref_id;
            }

            // Slice data blocks
            for (j = 0; j < num_blocks; j++) {
                // read and discard, unless it's the ref-ID block
                if (!(blk = cram_read_block(in_c)))
                    goto err;
                if (cram_block_get_content_id(blk) == embed_id) {
                    cram_uncompress_block(blk);
                    //printf("%.*s\n", blk->uncomp_size, blk->data);

                    int32_t usize = cram_block_get_uncomp_size(blk);
                    int ref_end = ref_start + usize;
                    if (ref_end > ref_len+1)
                        ref_end = ref_len+1;
                    if (ref_end > ref_start)
                        memcpy(ref + ref_start-1, cram_block_get_data(blk),
                               ref_end - ref_start);
                }
                cram_free_block(blk);
                blk = NULL;
            }
            cram_free_slice_header(shdr);
            shdr = NULL;
        }

        cram_free_container(c);
        c = NULL;
    }

    int ret = 0;
    if (curr_ref_id >= 0) {
        ret = dump_ref(h, iter, curr_ref_id, ref, ref_len, outfp, verbose);
    } else if (reg) {
        // no data present
        // no data present, but we explicitly asked for the reference so
        // report it still as Ns.
        ref_len = MIN(iter->end,  sam_hdr_tid2len(h, iter->tid));
        ref = malloc(ref_len);
        memset(ref, 'N', ref_len);
        if (!ref)
            goto err;
        ret = dump_ref(h, iter, iter->tid, ref, ref_len, outfp, verbose);
    }

    free(ref);
    if (iter)
        hts_itr_destroy(iter);

    return ret;

 err:
    free(ref);
    if (blk)
        cram_free_block(blk);
    if (shdr)
        cram_free_slice_header(shdr);
    if (c)
        cram_free_container(c);
    if (iter)
        hts_itr_destroy(iter);

    return -1;
}

/*
 * ---------------------------------------------------------------------------
 * MD method of reference construction
 */

// Returns the next cigar op code: one of the BAM_C* codes,
// or -1 if no more are present.
static inline
int next_cigar_op(uint32_t *cigar, int *ncigar, int *skip, int *spos,
                  uint32_t *cig_ind, uint32_t *cig_op, uint32_t *cig_len) {
    for(;;) {
        while (*cig_len == 0) {
            if (*cig_ind < *ncigar) {
                *cig_op  = cigar[*cig_ind] & BAM_CIGAR_MASK;
                *cig_len = cigar[*cig_ind] >> BAM_CIGAR_SHIFT;
                (*cig_ind)++;
            } else {
                return -1;
            }
        }

        if (skip[*cig_op]) {
            *spos += (bam_cigar_type(*cig_op)&1) * *cig_len;
            *cig_len = 0;
            continue;
        }

        (*cig_len)--;
        break;
    }

    return *cig_op;
}

// Converts a bam object with SEQ, POS/CIGAR and MD:Z to a reference.
// Updates ref[] array.
//
// Returns >0 on success,
//          0 on no-MD found,
//         -1 on failure (eg inconsistent data)
static int build_ref(bam1_t *b, char *ref, size_t ref_len) {
    uint8_t *seq = bam_get_seq(b);
    uint32_t *cigar = bam_get_cigar(b);
    int ncigar = b->core.n_cigar;
    uint32_t cig_op = 0, cig_len = 0, cig_ind = 0;

    const uint8_t *MD = bam_aux_get(b, "MD");
    if (!MD || *MD != 'Z')
        return 0;
    MD++;

    // Walk through MD + seq to generate ref
    int iseq = 0, iref = b->core.pos, next_op;
    int cig_skip[16] = {0,1,0,1,1,1,1,0,0,1,1,1,1,1,1,1};
    while (iseq < b->core.l_qseq && *MD) {
        if (isdigit(*MD)) {
            // match
            int len = strtol((char *)MD, (char **)&MD, 10);
            while (iseq < b->core.l_qseq && len) {
                if ((next_op = next_cigar_op(cigar, &ncigar, cig_skip,
                                             &iseq, &cig_ind, &cig_op,
                                             &cig_len)) < 0)
                    return -1;

                if (next_op != BAM_CMATCH &&
                    next_op != BAM_CEQUAL) {
                    print_error("MD2ref",
                                "MD:Z and CIGAR are incompatible");
                    return -1;
                }

                if (iref < ref_len)
                    ref[iref] = seq_nt16_str[bam_seqi(seq, iseq)];
                iseq++;
                iref++;
                len--;
            }
        } else if (*MD == '^') {
            // deletion
            MD++;
            while (*MD && isalpha(*MD)) {
                if ((next_op = next_cigar_op(cigar, &ncigar, cig_skip,
                                             &iseq, &cig_ind, &cig_op,
                                             &cig_len)) < 0)
                    return -1;

                if (next_op != BAM_CDEL) {
                    print_error("MD2ref",
                                "MD:Z and CIGAR are incompatible");
                    return -1;
                }

                if (iref < ref_len)
                    ref[iref] = *MD;

                MD++;
                iref++;
            }
        } else {
            // substitution
            if ((next_op = next_cigar_op(cigar, &ncigar, cig_skip,
                                         &iseq, &cig_ind, &cig_op,
                                         &cig_len)) < 0)
                return -1;

            if (next_op != BAM_CMATCH && next_op != BAM_CDIFF) {
                print_error("MD2ref", "MD:Z and CIGAR are incompatible");
                return -1;
            }
            if (iref < ref_len)
                ref[iref] = *MD;

            MD++;
            iref++;
            iseq++;
        }
    }

    return 1;
}

static int MD2ref(samFile *in, sam_hdr_t *h, hts_idx_t *idx, char *reg,
                  FILE *outfp, int verbose) {
    bam1_t *b = bam_init1();
    int r, last_tid = -99;
    size_t ref_len = 0;
    char *ref = NULL;
    int ret = -1;

    hts_itr_t *iter = NULL;
    if (idx && reg) {
        iter = sam_itr_querys(idx, h, reg);
        if (!iter) {
            print_error("reference", "failed to parse region '%s'", reg);
            goto err;
        }
    }

    while ((r = iter
                ? sam_itr_next(in, iter, b)
                : sam_read1(in, h, b)) >= 0) {
        // check b->core.tid and flush old seq.
        if (b->core.tid != last_tid) {
            if (last_tid >= 0)
                if (dump_ref(h, iter, last_tid, ref, ref_len, outfp,
                             verbose) < 0)
                    goto err;

            last_tid = b->core.tid;
            ref_len = sam_hdr_tid2len(h, last_tid);
            if (ref_len) {
                char *ref2 = realloc(ref, ref_len);
                if (!ref2)
                    goto err;
                else
                    ref = ref2;
                memset(ref, 'N', ref_len);
            }
        }

        if (build_ref(b, ref, ref_len) < 0)
            goto err;
    }

    if (last_tid >= 0) {
        if (dump_ref(h, iter, last_tid, ref, ref_len, outfp, verbose) < 0)
            goto err;
    } else if (reg) {
        // no data present, but we explicitly asked for the reference so
        // report it still as Ns.
        ref_len = MIN(iter->end,  sam_hdr_tid2len(h, iter->tid));
        ref = malloc(ref_len);
        memset(ref, 'N', ref_len);
        if (!ref)
            goto err;
        if (dump_ref(h, iter, iter->tid, ref, ref_len, outfp, verbose) < 0)
            goto err;
    }

    if (r < -1)
        goto err;

    ret = 0;

 err:
    if (iter)
        hts_itr_destroy(iter);
    bam_destroy1(b);
    free(ref);
    return ret;
}

int main_reference(int argc, char *argv[])
{
    int c, usage = 0, verbose = 1, use_embedded = 0;
    sam_hdr_t *h = 0;
    samFile *in = NULL;
    hts_idx_t *idx = NULL;
    sam_global_args ga;
    FILE *outfp = stdout;
    char *reg = NULL;

    static const struct option lopts[] = {
        {"output",   required_argument, NULL, 'o'},
        {"quiet",    no_argument,       NULL, 'q'},
        {"embedded", no_argument,       NULL, 'e'},
        {"region",   required_argument, NULL, 'r'},
        SAM_OPT_GLOBAL_OPTIONS('-', '-', '-', '-', '-', '@'),
        { NULL, 0, NULL, 0 }
    };

    sam_global_args_init(&ga);

    while ((c = getopt_long(argc, argv, "@:qo:er:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'o':
            if (!(outfp = fopen(optarg, "w"))) {
                perror(optarg);
                goto err;
            }
            break;

        case 'q':
            verbose = 0;
            break;

        case 'e':
            use_embedded = 1;
            break;

        case 'r':
            reg = optarg;
            break;

        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?': usage=1; break;
        }
    }

    if ((optind == argc && isatty(0)) || usage) {
        printf("Usage: samtools reference [-@ N] [-r region] [-e] [-q] [-o out.fa] [in.cram]\n");
        return 0;
    }

    char *fn = optind < argc ? argv[optind] : "-";
    if (!(in = sam_open(fn, "r"))) {
        print_error_errno("reference", "failed to open file '%s'", fn);
        return 1;
    }

    if (ga.nthreads > 0)
        hts_set_threads(in, ga.nthreads);

    if (!(h = sam_hdr_read(in)))
        goto err;

    if (reg) {
        idx = sam_index_load(in, fn);
        if (!idx) {
            print_error_errno("reference", "Failed to load the index");
            goto err;
        }
    }

    int ret = use_embedded
        ? cram2ref(in, h, idx, reg, outfp, verbose)
        : MD2ref(in, h, idx, reg, outfp, verbose);

    sam_hdr_destroy(h);
    if (outfp != stdout)
        fclose(outfp);
    if (idx)
        hts_idx_destroy(idx);
    sam_close(in);

    return ret;

 err:
    if (idx)
        hts_idx_destroy(idx);
    if (in)
        sam_close(in);
    if (h)
        sam_hdr_destroy(h);

    return 1;
}
