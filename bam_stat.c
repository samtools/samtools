/*  bam_stat.c -- flagstat subcommand.

    Copyright (C) 2009, 2011, 2013-2015, 2019 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <getopt.h>

#include "htslib/sam.h"
#include "samtools.h"
#include "sam_opts.h"

typedef struct {
    long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
    long long n_sgltn[2], n_read1[2], n_read2[2];
    long long n_dup[2];
    long long n_diffchr[2], n_diffhigh[2];
    long long n_secondary[2], n_supp[2];
} bam_flagstat_t;

#define flagstat_loop(s, c) do {                                        \
        int w = ((c)->flag & BAM_FQCFAIL)? 1 : 0;                       \
        ++(s)->n_reads[w];                                              \
        if ((c)->flag & BAM_FSECONDARY ) {                              \
            ++(s)->n_secondary[w];                                      \
        } else if ((c)->flag & BAM_FSUPPLEMENTARY ) {                   \
            ++(s)->n_supp[w];                                           \
        } else if ((c)->flag & BAM_FPAIRED) {                           \
            ++(s)->n_pair_all[w];                                       \
            if (((c)->flag & BAM_FPROPER_PAIR) && !((c)->flag & BAM_FUNMAP) ) ++(s)->n_pair_good[w];    \
            if ((c)->flag & BAM_FREAD1) ++(s)->n_read1[w];              \
            if ((c)->flag & BAM_FREAD2) ++(s)->n_read2[w];              \
            if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP)) ++(s)->n_sgltn[w];  \
            if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP)) { \
                ++(s)->n_pair_map[w];                                   \
                if ((c)->mtid != (c)->tid) {                            \
                    ++(s)->n_diffchr[w];                                \
                    if ((c)->qual >= 5) ++(s)->n_diffhigh[w];           \
                }                                                       \
            }                                                           \
        }                                                               \
        if (!((c)->flag & BAM_FUNMAP)) ++(s)->n_mapped[w];              \
        if ((c)->flag & BAM_FDUP) ++(s)->n_dup[w];                      \
    } while (0)

bam_flagstat_t *bam_flagstat_core(samFile *fp, sam_hdr_t *h)
{
    bam_flagstat_t *s;
    bam1_t *b;
    bam1_core_t *c;
    int ret;
    s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));
    b = bam_init1();
    c = &b->core;
    while ((ret = sam_read1(fp, h, b)) >= 0)
        flagstat_loop(s, c);
    bam_destroy1(b);
    if (ret != -1)
        fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");
    return s;
}

static const char *percent(char *buffer, long long n, long long total)
{
    if (total != 0) sprintf(buffer, "%.2f%%", (float)n / total * 100.0);
    else strcpy(buffer, "N/A");
    return buffer;
}

static const char *percent_json(char *buffer, long long n, long long total)
{
    if (total != 0) sprintf(buffer, "%.2f", (float)n / total * 100.0);
    else strcpy(buffer, "null");
    return buffer;
}

static void usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools flagstat [options] <in.bam>\n");
    sam_global_opt_help(fp, "-.---@-.");
    fprintf(fp, "  -O, --");
    fprintf(fp, "output-fmt FORMAT[,OPT[=VAL]]...\n"
            "               Specify output format (json, tsv)\n");
    exit(exit_status);
}

static void out_fmt_default(bam_flagstat_t *s)
{
    char b0[16], b1[16];
    printf("%lld + %lld in total (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
    printf("%lld + %lld secondary\n", s->n_secondary[0], s->n_secondary[1]);
    printf("%lld + %lld supplementary\n", s->n_supp[0], s->n_supp[1]);
    printf("%lld + %lld duplicates\n", s->n_dup[0], s->n_dup[1]);
    printf("%lld + %lld mapped (%s : %s)\n", s->n_mapped[0], s->n_mapped[1], percent(b0, s->n_mapped[0], s->n_reads[0]), percent(b1, s->n_mapped[1], s->n_reads[1]));
    printf("%lld + %lld paired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
    printf("%lld + %lld read1\n", s->n_read1[0], s->n_read1[1]);
    printf("%lld + %lld read2\n", s->n_read2[0], s->n_read2[1]);
    printf("%lld + %lld properly paired (%s : %s)\n", s->n_pair_good[0], s->n_pair_good[1], percent(b0, s->n_pair_good[0], s->n_pair_all[0]), percent(b1, s->n_pair_good[1], s->n_pair_all[1]));
    printf("%lld + %lld with itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
    printf("%lld + %lld singletons (%s : %s)\n", s->n_sgltn[0], s->n_sgltn[1], percent(b0, s->n_sgltn[0], s->n_pair_all[0]), percent(b1, s->n_sgltn[1], s->n_pair_all[1]));
    printf("%lld + %lld with mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
    printf("%lld + %lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
}

static void out_fmt_json(bam_flagstat_t *s) {
    char b0[16], b1[16];
    printf("{\n \"QC-passed reads\": { \n"
                 "  \"total\": %lld, \n"
                 "  \"secondary\": %lld, \n"
                 "  \"supplementary\": %lld, \n"
                 "  \"duplicates\": %lld, \n"
                 "  \"mapped\": %lld, \n"
                 "  \"mapped %%\": %s, \n"
                 "  \"paired in sequencing\": %lld, \n"
                 "  \"read1\": %lld, \n"
                 "  \"read2\": %lld, \n"
                 "  \"properly paired\": %lld, \n"
                 "  \"properly paired %%\": %s, \n"
                 "  \"with itself and mate mapped\": %lld, \n"
                 "  \"singletons\": %lld, \n"
                 "  \"singletons %%\": %s, \n"
                 "  \"with mate mapped to a different chr\": %lld, \n"
                 "  \"with mate mapped to a different chr (mapQ >= 5)\": %lld \n"
                 " },"
            "\n \"QC-failed reads\": { \n"
                 "  \"total\": %lld, \n"
                 "  \"secondary\": %lld, \n"
                 "  \"supplementary\": %lld, \n"
                 "  \"duplicates\": %lld, \n"
                 "  \"mapped\": %lld, \n"
                 "  \"mapped %%\": %s, \n"
                 "  \"paired in sequencing\": %lld, \n"
                 "  \"read1\": %lld, \n"
                 "  \"read2\": %lld, \n"
                 "  \"properly paired\": %lld, \n"
                 "  \"properly paired %%\": %s, \n"
                 "  \"with itself and mate mapped\": %lld, \n"
                 "  \"singletons\": %lld, \n"
                 "  \"singletons %%\": %s, \n"
                 "  \"with mate mapped to a different chr\": %lld, \n"
                 "  \"with mate mapped to a different chr (mapQ >= 5)\": %lld \n"
                 " }\n"
            "}\n",
        s->n_reads[0],
        s->n_secondary[0],
        s->n_supp[0],
        s->n_dup[0],
        s->n_mapped[0],
        percent_json(b0, s->n_mapped[0], s->n_reads[0]),
        s->n_pair_all[0],
        s->n_read1[0],
        s->n_read2[0],
        s->n_pair_good[0],
        percent_json(b0, s->n_pair_good[0], s->n_pair_all[0]),
        s->n_pair_map[0],
        s->n_sgltn[0],
        percent_json(b0, s->n_sgltn[0], s->n_pair_all[0]),
        s->n_diffchr[0],
        s->n_diffhigh[0],
        s->n_reads[1],
        s->n_secondary[1],
        s->n_supp[1],
        s->n_dup[1],
        s->n_mapped[1],
        percent_json(b1, s->n_mapped[1], s->n_reads[1]),
        s->n_pair_all[1],
        s->n_read1[1],
        s->n_read2[1],
        s->n_pair_good[1],
        percent_json(b1, s->n_pair_good[1], s->n_pair_all[1]),
        s->n_pair_map[1],
        s->n_sgltn[1],
        percent_json(b1, s->n_sgltn[1], s->n_pair_all[1]),
        s->n_diffchr[1],
        s->n_diffhigh[1]
    );
}

static void out_fmt_tsv(bam_flagstat_t *s) {
    char b0[16], b1[16];
    printf("%lld\t%lld\ttotal (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
    printf("%lld\t%lld\tsecondary\n", s->n_secondary[0], s->n_secondary[1]);
    printf("%lld\t%lld\tsupplementary\n", s->n_supp[0], s->n_supp[1]);
    printf("%lld\t%lld\tduplicates\n", s->n_dup[0], s->n_dup[1]);
    printf("%lld\t%lld\tmapped\n", s->n_mapped[0], s->n_mapped[1]);
    printf("%s\t%s\tmapped %%\n", percent(b0, s->n_mapped[0], s->n_reads[0]), percent(b1, s->n_mapped[1], s->n_reads[1]));
    printf("%lld\t%lld\tpaired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
    printf("%lld\t%lld\tread1\n", s->n_read1[0], s->n_read1[1]);
    printf("%lld\t%lld\tread2\n", s->n_read2[0], s->n_read2[1]);
    printf("%lld\t%lld\tproperly paired\n", s->n_pair_good[0], s->n_pair_good[1]);
    printf("%s\t%s\tproperly paired %%\n", percent(b0, s->n_pair_good[0], s->n_pair_all[0]), percent(b1, s->n_pair_good[1], s->n_pair_all[1]));
    printf("%lld\t%lld\twith itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
    printf("%lld\t%lld\tsingletons\n", s->n_sgltn[0], s->n_sgltn[1]);
    printf("%s\t%s\tsingletons %%\n", percent(b0, s->n_sgltn[0], s->n_pair_all[0]), percent(b1, s->n_sgltn[1], s->n_pair_all[1]));
    printf("%lld\t%lld\twith mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
    printf("%lld\t%lld\twith mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
}

/*
 * Select flagstats output format to print.
 */
static void output_fmt(bam_flagstat_t *s, const char *out_fmt)
{
  if (strcmp(out_fmt, "json") == 0 || strcmp(out_fmt, "JSON") == 0) {
    out_fmt_json(s);
  } else if (strcmp(out_fmt, "tsv") == 0 || strcmp(out_fmt, "TSV") == 0) {
    out_fmt_tsv(s);
  } else {
    out_fmt_default(s);
  }
}

int bam_flagstat(int argc, char *argv[])
{
    samFile *fp;
    sam_hdr_t *header;
    bam_flagstat_t *s;
    const char *out_fmt = "default";
    int c;

    enum {
        INPUT_FMT_OPTION = CHAR_MAX+1,
    };

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', '-', '-', '@'),
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "@:O:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'O':
          out_fmt = optarg;
          break;
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
        print_error_errno("flagstat", "Cannot open input file \"%s\"", argv[optind]);
        return 1;
    }
    if (ga.nthreads > 0)
        hts_set_threads(fp, ga.nthreads);

    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS,
                    SAM_FLAG | SAM_MAPQ | SAM_RNEXT)) {
        fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
        return 1;
    }

    if (hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        return 1;
    }

    header = sam_hdr_read(fp);
    if (header == NULL) {
        fprintf(stderr, "Failed to read header for \"%s\"\n", argv[optind]);
        return 1;
    }

    s = bam_flagstat_core(fp, header);
    output_fmt(s, out_fmt);
    free(s);
    sam_hdr_destroy(header);
    sam_close(fp);
    sam_global_args_free(&ga);
    return 0;
}
