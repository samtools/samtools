/*  bam_checksum.c -- produces checksums on SAM/BAM/CRAM/FASTA/FASTQ data

    Copyright (C) 2024 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

The primary work here is GRL since 2021, under an MIT license.
Sections derived from Gap5, which include calculate_consensus_gap5()
associated functions, are mostly copyright Genome Research Limited from
2003 onwards.  These were originally under a BSD license, but as GRL is
copyright holder these portions can be considered to also be under the
same MIT license below:


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

/*
 * This is inspired by Biobambam's bamseqchksum tool and ideas from
 * German Tischler and David Jackson.
 *
 * It computes order agnostic checksums for a variety of SAM fields, allowing
 * validation that all the data is still present at different stages of an
 * analysis pipeline.  This may be useful to detect sequences which have been
 * lost by an aligner, memory corruptions flipping individual sequence bases,
 * or file format decoding errors.
 *
 * We start with something basic such as a FASTQ file, and name, seq and qual
 * checksums should still all match after aligning and sorting.
 */

/*
TODO

- Support multiple read-groups, which aids spliting pooled samples and
  tracking that data isn't lost.
- Separate "all" from "pass" only (dropping QC fail)
- Mechanisms for merging checksums together.  Eg merging two read-groups into
  a new file.
- Make tags configurable
 */

#include <config.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>

#include <htslib/sam.h>
#include "sam_opts.h"
#include "sam_utils.h"

// TODO: consider adding hts_crc32 to htslib, similar to our hts_md5
// This exposes libdeflate from htslib instead of punting the configure /
// building requirements downstream.
#ifdef HAVE_LIBDEFLATE
#  include <libdeflate.h>
#  define crc32 libdeflate_crc32
#endif


typedef struct {
    int incl_flags, req_flags, excl_flags;  // BAM flags filtering
    int nthreads;
    int rev_comp;
} opts;

/*
 * The hash is multiplicative within a finite field, modulo PRIME.
 * We need to avoid zeros, and the data type has to be large enough to ensure
 * no wraparound happens (other than the intended modulo).
 */
#define PRIME ((1u<<31)-1)
#ifdef HASH_ADD
uint64_t update_hash(uint64_t hash, uint32_t crc) {
    return (hash + crc) % PRIME;
}
#else
uint64_t update_hash(uint64_t hash, uint32_t crc) {
    crc &= PRIME;
    if (crc == 0 || crc == PRIME)
	crc = 1;

    return (hash * crc) % PRIME;
}
#endif

int checksum(sam_global_args *ga, opts *o, char *fn) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    uint8_t *seq_buf = NULL;
    uint8_t *qual_buf = NULL;
    bam1_t *b = bam_init1();
    static const char *tags[] = {"BC","FI","QT","RT","TC"};
    const int ntags = sizeof(tags) / sizeof(*tags);
    uint8_t **tag_ptr = calloc(ntags, sizeof(*tag_ptr));
    size_t   *tag_len = calloc(ntags, sizeof(*tag_len));

    if (!b || !tag_ptr || !tag_len)
	goto err;

    int tag_keep[125][125] = {0};
    for (int i = 0; i < ntags; i++)
	tag_keep[tags[i][0]-'0'][tags[i][1]-'0'] = i+1;

#ifdef HASH_ADD
    uint64_t name_hash = 0;
    uint64_t seq_hash  = 0;
    uint64_t qual_hash = 0;
    uint64_t aux_hash  = 0;
#else
    uint64_t name_hash = 1;
    uint64_t seq_hash  = 1;
    uint64_t qual_hash = 1;
    uint64_t aux_hash  = 1;
#endif

    const uint32_t crc32_start = crc32(0L, Z_NULL, 0);

    fp = sam_open_format(fn, "r", &ga->in);
    if (!fp) {
	print_error_errno("checksum", "Cannot open input file \"%s\"", fn);
	goto err;
    }

    if (!(hdr = sam_hdr_read(fp)))
	goto err;

    // FIXME: use kstring instead
    uint8_t *aux_dat = NULL;
    size_t aux_sz = 0;
    size_t seq_buf_len = 0;
    int r;
    uint64_t count = 0;

    while ((r = sam_read1(fp, hdr, b)) >= 0) {
	// TODO: configurable filter
	if (b->core.flag & o->excl_flags)
	    continue;

	// 8 bits of flag corresponding to original instrument data
	uint8_t flags = b->core.flag & (BAM_FPAIRED | BAM_FREAD1 | BAM_FREAD2);

	// Copy sequence out from nibble to base, and reverse complement
	// seq / qual if required.  Qual is +33 (ASCII format) only for
	// compatibility with biobambam's bamseqchksum tool.
	uint8_t *seq = bam_get_seq(b);
	uint8_t *qual = bam_get_qual(b);
	if (seq_buf_len < b->core.l_qseq) {
	    uint8_t *tmp;
	    
	    if (!(tmp = realloc(seq_buf, seq_buf_len = b->core.l_qseq)))
		goto err;
	    seq_buf = tmp;

	    if (!(tmp = realloc(qual_buf, seq_buf_len)))
		goto err;
	    qual_buf = tmp;
	}
	if ((b->core.flag & BAM_FREVERSE) && o->rev_comp) {
	    for (int i=0, j=b->core.l_qseq-1; i < b->core.l_qseq; i++,j--) {
		seq_buf[j] = "=TGKCYSBAWRDMHVN"[bam_seqi(seq, i)];
		qual_buf[j] = qual[i]+33;
	    }
	    qual = qual_buf;
	} else {
	    for (int i = 0; i < b->core.l_qseq; i++) {
		seq_buf[i] = seq_nt16_str[bam_seqi(seq, i)];
		qual_buf[i] = qual[i]+33;
	    }
	    qual = qual_buf;
	}

	// name + flag + seq.
	// Name includes single nul byte, for compatibility with bamseqchksum.
	// (Would be better to be flag + seq + name as they wouldn't need to
	// recompute CRC for flag+seq twice.)
	uint32_t crc = crc32(crc32_start, (uint8_t *)bam_get_qname(b),
			     b->core.l_qname - b->core.l_extranul);
	crc = crc32(crc, &flags, 1);
	crc = crc32(crc, seq_buf, b->core.l_qseq);
	name_hash = update_hash(name_hash, crc);

	// flag + seq
	crc = crc32(crc32_start, &flags, 1);
	crc = crc32(crc, seq_buf, b->core.l_qseq);
	seq_hash = update_hash(seq_hash, crc);
	uint32_t crc_flag_seq = crc;

	// flag + seq + qual
	crc = crc32(crc_flag_seq, qual, b->core.l_qseq);
	qual_hash = update_hash(qual_hash, crc); 

	// flag + seq + aux tags
	size_t aux_len = bam_get_l_aux(b);
	if (aux_sz < aux_len)
	    aux_dat = realloc(aux_dat, aux_sz = aux_len); // FIXME kstring
	uint8_t *aux_ptr = aux_dat;

	// Pass 1: find all tags to copy and their lengths
	uint8_t *aux = bam_aux_first(b), *aux_next;
	memset(tag_len, 0, ntags * sizeof(*tag_len));
	while (aux) {
	    aux_next = bam_aux_next(b, aux);
	    if (!(aux[-2] >= '0' && aux[-2] <= 'z' &&
		  aux[-1] >= '0' && aux[-2] <= 'z'))
		continue; // skip illegal tag names
	    int i = tag_keep[aux[-2]-'0'][aux[-1]-'0']-1;
	    if (i>=0) {
		// found one
		size_t tag_sz = aux_next
		    ? aux_next - aux
		    : b->data + b->l_data - aux + 2;

		tag_ptr[i] = aux-2;
		tag_len[i] = tag_sz;
	    }

	    aux = aux_next;
	}

	// Pass 2: copy tags in the order we requested
	for (int i = 0; i < ntags; i++) {
	    if (tag_len[i]) {
		memcpy(aux_ptr, tag_ptr[i], tag_len[i]);
		aux_ptr += tag_len[i];
	    }
	}

	if (aux_ptr > aux_dat) {
	    crc = crc32(crc_flag_seq, aux_dat, aux_ptr - aux_dat);
	    aux_hash = update_hash(aux_hash, crc);
	} else {
	    aux_hash = update_hash(aux_hash, crc_flag_seq);
	}

	count++;
    }

    printf("Count          %ld\n", count);
    printf("Flag+Seq       %08lx\n", seq_hash);
    printf("Name+Flag+Seq  %08lx\n", name_hash);
    printf("Flag+Seq+Qual  %08lx\n", qual_hash);
    printf("Flag+Seq+Aux   %08lx\n", aux_hash);
    puts("");

    if (r <= -1)
	goto err;
    if (hdr)
	sam_hdr_destroy(hdr);

    if (sam_close(fp) < 0) {
	print_error_errno("checksum", "Closing input file \"%s\"", fn);
	goto err;
    }

    free(seq_buf);
    free(qual_buf);
    free(tag_ptr);
    free(tag_len);

    bam_destroy1(b);
    return 0;

 err:
    if (b)   bam_destroy1(b);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp)  sam_close(fp);

    free(seq_buf);
    free(qual_buf);
    free(tag_ptr);
    free(tag_len);

    return -1;
}

void usage_exit(FILE *fp, int ret) {
    fprintf(stderr, "Usage: samtools checksum [options] [file]\n");
    exit(ret);
}

int main_checksum(int argc, char **argv) {
    opts opts = {
	.incl_flags   = 0xffff,
	.req_flags    = 0,
	.excl_flags   = BAM_FSECONDARY | BAM_FSUPPLEMENTARY,
	.rev_comp     = 1,
    };

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 'I', '-', '-', '.', '@'),
	{"--excl-flags",    required_argument, NULL, 'F'},
	{"--exclude-flags", required_argument, NULL, 'F'},
	{"--require-flags", required_argument, NULL, 'f'},
	{"--incl-flags",    required_argument, NULL, 1},
	{"--include-flags", required_argument, NULL, 1},
	{"--threads",       required_argument, NULL, '@'},
	{NULL, 0, NULL, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "@:f:F:", lopts, NULL)) >= 0) {
	switch (c) {
	case '@': opts.nthreads = atoi(optarg); break;
	case 'F': opts.excl_flags = atoi(optarg); break;
	case 'f': opts.req_flags = atoi(optarg); break;
	case  1 : opts.incl_flags = atoi(optarg); break;
	default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
	    /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
	}
    }

    int ret = 0;
    if (argc-optind) {
	while (optind < argc)
	    ret |= checksum(&ga, &opts, argv[optind++]) < 0;
    } else {
	ret = checksum(&ga, &opts, "-") < 0;
    }

    return ret;
}
