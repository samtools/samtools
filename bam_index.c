/*  bam_index.c -- index and idxstats subcommands.

    Copyright (C) 2008-2011, 2013-2016, 2018, 2019  Genome Research Ltd.
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

#include <config.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <stdlib.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>
#include <getopt.h>

#include "samtools.h"
#include "sam_opts.h"

#define BAM_LIDX_SHIFT    14

static void index_usage(FILE *fp)
{
    fprintf(fp,
"Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]\n"
"Options:\n"
"  -b       Generate BAI-format index for BAM files [default]\n"
"  -c       Generate CSI-format index for BAM files\n"
"  -m INT   Set minimum interval size for CSI indices to 2^INT [%d]\n"
"  -@ INT   Sets the number of threads [none]\n", BAM_LIDX_SHIFT);
}

int bam_index(int argc, char *argv[])
{
    int csi = 0;
    int min_shift = BAM_LIDX_SHIFT;
    int n_threads = 0;
    int c, ret;

    while ((c = getopt(argc, argv, "bcm:@:")) >= 0)
        switch (c) {
        case 'b': csi = 0; break;
        case 'c': csi = 1; break;
        case 'm': csi = 1; min_shift = atoi(optarg); break;
        case '@': n_threads = atoi(optarg); break;
        default:
            index_usage(stderr);
            return 1;
        }

    if (optind == argc) {
        index_usage(stdout);
        return 1;
    }

    ret = sam_index_build3(argv[optind], argv[optind+1], csi? min_shift : 0, n_threads);
    switch (ret) {
    case 0:
        return 0;

    case -2:
        print_error_errno("index", "failed to open \"%s\"", argv[optind]);
        break;

    case -3:
        print_error("index", "\"%s\" is in a format that cannot be usefully indexed", argv[optind]);
        break;

    case -4:
        if (argv[optind+1])
            print_error("index", "failed to create or write index \"%s\"", argv[optind+1]);
        else
            print_error("index", "failed to create or write index");
        break;

    default:
        print_error_errno("index", "failed to create index for \"%s\"", argv[optind]);
        break;
    }

    return EXIT_FAILURE;
}

/*
 * Cram indices do not contain mapped/unmapped record counts, so we have to
 * decode each record and count.  However we can speed this up as much as
 * possible by using the required fields parameter.
 *
 * This prints the stats to stdout in the same manner than the BAM function
 * does.
 *
 * Returns 0 on success,
 *        -1 on failure.
 */
int slow_idxstats(samFile *fp, sam_hdr_t *header) {
    int ret, last_tid = -2;
    bam1_t *b = bam_init1();

    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, SAM_RNAME | SAM_FLAG))
        return -1;

    uint64_t (*count0)[2] = calloc(sam_hdr_nref(header)+1, sizeof(*count0));
    uint64_t (*counts)[2] = count0+1;
    if (!count0)
        return -1;

    while ((ret = sam_read1(fp, header, b)) >= 0) {
        if (b->core.tid >= sam_hdr_nref(header) || b->core.tid < -1) {
            free(count0);
            return -1;
        }

        if (b->core.tid != last_tid) {
            if (last_tid >= -1) {
                if (counts[b->core.tid][0] + counts[b->core.tid][1]) {
                    print_error("idxstats", "file is not position sorted");
                    free(count0);
                    return -1;
                }
            }
            last_tid = b->core.tid;
        }

        counts[b->core.tid][(b->core.flag & BAM_FUNMAP) ? 1 : 0]++;
    }

    if (ret == -1) {
        int i;
        for (i = 0; i < sam_hdr_nref(header); i++) {
            printf("%s\t%"PRId64"\t%"PRIu64"\t%"PRIu64"\n",
                   sam_hdr_tid2name(header, i),
                   (int64_t) sam_hdr_tid2len(header, i),
                   counts[i][0], counts[i][1]);
        }
        printf("*\t0\t%"PRIu64"\t%"PRIu64"\n", counts[-1][0], counts[-1][1]);
    }

    free(count0);

    bam_destroy1(b);

    return (ret == -1) ? 0 : -1;
}

static void usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools idxstats [options] <in.bam>\n");
    sam_global_opt_help(fp, "-.---@-.");
    exit(exit_status);
}

int bam_idxstats(int argc, char *argv[])
{
    hts_idx_t* idx;
    sam_hdr_t* header;
    samFile* fp;
    int c;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', '-', '@'),
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "@:", lopts, NULL)) >= 0) {
        switch (c) {
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    if (argc != optind+1) {
        if (argc == optind) usage_exit(stdout, EXIT_SUCCESS);
        else usage_exit(stderr, EXIT_FAILURE);
    }

    fp = sam_open_format(argv[optind], "r", &ga.in);
    if (fp == NULL) {
        print_error_errno("idxstats", "failed to open \"%s\"", argv[optind]);
        return 1;
    }
    header = sam_hdr_read(fp);
    if (header == NULL) {
        print_error("idxstats", "failed to read header for \"%s\"", argv[optind]);
        return 1;
    }

    if (hts_get_format(fp)->format != bam) {
    slow_method:
        if (ga.nthreads)
            hts_set_threads(fp, ga.nthreads);

        if (slow_idxstats(fp, header) < 0) {
            print_error("idxstats", "failed to process \"%s\"", argv[optind]);
            return 1;
        }
    } else {
        idx = sam_index_load(fp, argv[optind]);
        if (idx == NULL) {
            print_error("idxstats", "fail to load index for \"%s\", "
                        "reverting to slow method", argv[optind]);
            goto slow_method;
        }

        int i;
        for (i = 0; i < sam_hdr_nref(header); ++i) {
            // Print out contig name and length
            printf("%s\t%"PRId64, sam_hdr_tid2name(header, i), (int64_t) sam_hdr_tid2len(header, i));
            // Now fetch info about it from the meta bin
            uint64_t u, v;
            hts_idx_get_stat(idx, i, &u, &v);
            printf("\t%" PRIu64 "\t%" PRIu64 "\n", u, v);
        }
        // Dump information about unmapped reads
        printf("*\t0\t0\t%" PRIu64 "\n", hts_idx_get_n_no_coor(idx));
        hts_idx_destroy(idx);
    }

    sam_hdr_destroy(header);
    sam_close(fp);
    return 0;
}
