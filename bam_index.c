/*  bam_index.c -- index and idxstats subcommands.

    Copyright (C) 2008-2011, 2013, 2014 Genome Research Ltd.
    Portions copyright (C) 2010 Broad Institute.
    Portions copyright (C) 2013 Peter Cock, The James Hutton Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <stdlib.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>

#define BAM_LIDX_SHIFT    14

int bam_index_build2(const char *fn, const char *_fnidx)
{
    fprintf(stderr, "Samtools-htslib-API: bam_index_build2() not yet implemented\n");
    abort();
}

static void index_usage(FILE *fp)
{
    fprintf(fp,
"Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]\n"
"Options:\n"
"  -b       Generate BAI-format index for BAM files [default]\n"
"  -c       Generate CSI-format index for BAM files\n"
"  -m INT   Set minimum interval size for CSI indices to 2^INT [%d]\n", BAM_LIDX_SHIFT);
}

int bam_index(int argc, char *argv[])
{
    int csi = 0;
    int min_shift = BAM_LIDX_SHIFT;
    int c;

    while ((c = getopt(argc, argv, "bcm:")) >= 0)
        switch (c) {
        case 'b': csi = 0; break;
        case 'c': csi = 1; break;
        case 'm': csi = 1; min_shift = atoi(optarg); break;
        default:
            index_usage(stderr);
            return 1;
        }

    if (optind == argc) {
        index_usage(stdout);
        return 1;
    }
    if (argc - optind > 1) bam_index_build2(argv[optind], argv[optind+1]);
    else bam_index_build(argv[optind], csi? min_shift : 0);
    return 0;
}

int bam_idxstats(int argc, char *argv[])
{
    hts_idx_t* idx;
    bam_hdr_t* header;
    samFile* fp;

    if (argc < 2) {
        fprintf(stderr, "Usage: samtools idxstats <in.bam>\n");
        return 1;
    }
    fp = sam_open(argv[1], "r");
    if (fp == NULL) { fprintf(stderr, "[%s] fail to open BAM.\n", __func__); return 1; }
    header = sam_hdr_read(fp);
    idx = sam_index_load(fp, argv[1]);
    if (idx == NULL) { fprintf(stderr, "[%s] fail to load the index.\n", __func__); return 1; }

    int i;
    for (i = 0; i < header->n_targets; ++i) {
        // Print out contig name and length
        printf("%s\t%d", header->target_name[i], header->target_len[i]);
        // Now fetch info about it from the meta bin
        uint64_t u, v;
        hts_idx_get_stat(idx, i, &u, &v);
        printf("\t%" PRIu64 "\t%" PRIu64 "\n", u, v);
    }
    // Dump information about unmapped reads
    printf("*\t0\t0\t%" PRIu64 "\n", hts_idx_get_n_no_coor(idx));
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    sam_close(fp);
    return 0;
}
