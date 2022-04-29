/*  cram_size.c -- produces summary of the size of each cram data-series

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
#include <strings.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/cram.h"
#include "htslib/kstring.h"
#include "samtools.h"
#include "sam_opts.h"

/*--------------------------------------------------
 * To add to htslib
 */
extern cram_block_slice_hdr *cram_decode_slice_header(cram_fd *fd, cram_block *b);
extern void cram_free_slice_header(cram_block_slice_hdr *hdr);

struct cram_block_slice_hdr {
    //enum cram_content_type content_type;
    int content_type;
    int32_t ref_seq_id;     /* if content_type == MAPPED_SLICE */
    int64_t ref_seq_start;  /* if content_type == MAPPED_SLICE */
    int64_t ref_seq_span;   /* if content_type == MAPPED_SLICE */
    int32_t num_records;
    int64_t record_counter;
    int32_t num_blocks;
    int32_t num_content_ids;
    int32_t *block_content_ids;
    int32_t ref_base_id;    /* if content_type == MAPPED_SLICE */
    unsigned char md5[16];
};

extern int cram_slice_hdr_get_num_blocks(cram_block_slice_hdr *h);

extern
cram_block_compression_hdr *cram_decode_compression_header(cram_fd *fd,
							   cram_block *b);
extern
void cram_free_compression_header(cram_block_compression_hdr *hdr);

/* External IDs used by this implementation (only assumed during writing) */
enum cram_DS_ID {
    DS_CORE   = 0,
    DS_aux    = 1, // aux_blk
    DS_aux_OQ = 2,
    DS_aux_BQ = 3,
    DS_aux_BD = 4,
    DS_aux_BI = 5,
    DS_aux_FZ = 6, // also ZM:B
    DS_aux_oq = 7, // other qualities
    DS_aux_os = 8, // other sequences
    DS_aux_oz = 9, // other strings
    DS_ref,
    DS_RN, // name_blk
    DS_QS, // qual_blk
    DS_IN, // base_blk
    DS_SC, // soft_blk

    DS_BF, // start loop
    DS_CF,
    DS_AP,
    DS_RG,
    DS_MQ,
    DS_NS,
    DS_MF,
    DS_TS,
    DS_NP,
    DS_NF,
    DS_RL,
    DS_FN,
    DS_FC,
    DS_FP,
    DS_DL,
    DS_BA,
    DS_BS,
    DS_TL,
    DS_RI,
    DS_RS,
    DS_PD,
    DS_HC,
    DS_BB,
    DS_QQ,

    DS_TN, // end loop

    DS_RN_len,
    DS_SC_len,
    DS_BB_len,
    DS_QQ_len,

    DS_TC, // CRAM v1.0 tags
    DS_TM, // test
    DS_TV, // test

    DS_END,
};

#define CRAM_MAP_HASH 32
#define CRAM_MAP(a,b) (((a)*3+(b))&(CRAM_MAP_HASH-1))
struct cram_block_compression_hdr {
    int32_t ref_seq_id;
    int64_t ref_seq_start;
    int64_t ref_seq_span;
    int32_t num_records;
    int32_t num_landmarks;
    int32_t *landmark;

    /* Flags from preservation map */
    int read_names_included;
    int AP_delta;
    // indexed by ref-base and subst. code
    char substitution_matrix[5][4];
    int no_ref;
    int qs_seq_orient; // 1 => same as seq. 0 => original orientation

    // TD Dictionary as a concatenated block
    cram_block *TD_blk;          // Tag Dictionary
    int nTL;                     // number of TL entries in TD
    unsigned char **TL;          // array of size nTL, pointer into TD_blk.
    //khash_t(m_s2i) *TD_hash;     // Keyed on TD strings, map to TL[] indices
    void *TD_hash;
    //string_alloc_t *TD_keys;     // Pooled keys for TD hash.
    void *TD_keys;

    //khash_t(map) *preservation_map;
    void *preservation_map;
    struct cram_map *rec_encoding_map[CRAM_MAP_HASH];
    struct cram_map *tag_encoding_map[CRAM_MAP_HASH];

    struct cram_codec *codecs[DS_END];

    char *uncomp; // A single block of uncompressed data
    size_t uncomp_size, uncomp_alloc;

    // Total codec count, used for index to block_by_id for transforms
    int ncodecs;
};

//struct cram_codec *cram_codec_for_DS(cram_block_compression_hdr *hdr,
//				     int DS, int aux) {
//    return hdr->comp_hdr->codecs[DS];
//}

/*
 * End of htslib additions
 *--------------------------------------------------
 */


/*=============================================================================
 * cram_size
 */

static int cram_size(samFile *in, sam_hdr_t *h, FILE *outfp, int verbose) {
    cram_fd *in_c;
    cram_container *c = NULL;
    cram_block *blk = NULL;
    cram_block_slice_hdr *shdr = NULL;

    in_c = in->fp.cram; // low level htslib abuse?
    while ((c = cram_read_container(in_c))) {
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

	// Container compression header
	int32_t num_slices;
	if (!(blk = cram_read_block(in_c)))
	    goto err;

	// Decode compression header...
	cram_block_compression_hdr *chdr;
	chdr = cram_decode_compression_header(in_c, blk);
	cram_free_block(blk);
	blk = NULL;

	cram_free_compression_header(chdr);
	
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

	    // Slice data blocks
	    for (j = 0; j < num_blocks; j++) {
		// read and discard, unless it's the ref-ID block
		if (!(blk = cram_read_block(in_c)))
		    goto err;

		int32_t csize = cram_block_get_comp_size(blk);
		int32_t usize = cram_block_get_uncomp_size(blk);
		int cid = cram_block_get_content_id(blk);

		fprintf(stderr, "Block id %d  csz %d  usz %d\n",
			cid, csize, usize);
		cram_free_block(blk);
		blk = NULL;
	    }
	    cram_free_slice_header(shdr);
	    shdr = NULL;
	}

	cram_free_container(c);
	c = NULL;
    }

    return 0;

 err:
    if (blk)
	cram_free_block(blk);
    if (shdr)
	cram_free_slice_header(shdr);
    if (c)
	cram_free_container(c);
    return -1;
}

int main_cram_size(int argc, char *argv[])
{
    int c, usage = 0, verbose = 1;
    sam_hdr_t *h = 0;
    samFile *in = NULL;
    sam_global_args ga;
    FILE *outfp = stdout;

    static const struct option lopts[] = {
        {"output", required_argument, NULL, 'o'},
        {"quiet",  no_argument, NULL, 'q'},
        SAM_OPT_GLOBAL_OPTIONS('-', '-', '-', '-', '-', '-'),
        { NULL, 0, NULL, 0 }
    };

    sam_global_args_init(&ga);

    while ((c = getopt_long(argc, argv, "qo:", lopts, NULL)) >= 0) {
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
	    
	default:
	    if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
	    /* else fall-through */
	case '?': usage=1; break;
	}
    }

    if ((optind == argc && isatty(0)) || usage) {
	printf("Usage: samtools cram_size [-o out.fa] [in.cram]\n");
	return 0;
    }

    char *fn = optind < argc ? argv[optind] : "-";

    if (!(in = sam_open(fn, "r"))) {
        print_error_errno("cram_size", "failed to open file '%s'", fn);
	return 1;
    }
    if (!(h = sam_hdr_read(in)))
	goto err;

    int ret = cram_size(in, h, outfp, verbose);
    sam_hdr_destroy(h);
    sam_close(in);
    if (outfp != stdout)
	fclose(outfp);

    return ret;

 err:
    if (in)
	sam_close(in);
    if (h)
	sam_hdr_destroy(h);
    
    return 1;
}
