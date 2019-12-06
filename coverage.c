/* coverage.c -- samtools coverage subcommand

    Copyright (C) 2018,2019 Florian Breitwieser
    Portions copyright (C) 2019 Genome Research Ltd.

    Author: Florian P Breitwieser <florian.bw@gmail.com>

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

/* This program calculates coverage from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bamcov -D_MAIN_BAMCOV coverage.c -lhts -lz
 */

// C headers
#include <config.h>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>  // variadic functions
#include <limits.h>  // INT_MAX
#include <math.h>    // round
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/ioctl.h>
#endif

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "samtools.h"
#include "sam_opts.h"

const char *VERSION = "0.1";

typedef struct {  // auxiliary data structure to hold a BAM file
    samFile *fp;     // file handle
    sam_hdr_t *hdr;  // file header
    hts_itr_t *iter; // iterator to a region - NULL for us by default
    int min_mapQ;    // mapQ filter
    int min_len;     // length filter
    unsigned int n_reads;  // records the number of reads seen in file
    unsigned int n_selected_reads; // records the number of reads passing filter
    unsigned long summed_mapQ; // summed mapQ of all reads passing filter
    int fail_flags;
    int required_flags;
} bam_aux_t;

typedef struct {  // auxiliary data structure to hold stats on coverage
    unsigned long long n_covered_bases;
    unsigned long long summed_coverage;
    unsigned long long summed_baseQ;
    unsigned long long summed_mapQ;
    unsigned int n_reads;
    unsigned int n_selected_reads;
    int32_t tid;    // chromosome ID, defined by header
    hts_pos_t beg;
    hts_pos_t end;
    int64_t bin_width;
} stats_aux_t;

#if __STDC_VERSION__ >= 199901L
#define VERTICAL_LINE "\u2502" // BOX DRAWINGS LIGHT VERTICAL

// UTF8 specifies block characters in eights going from \u2581 (lower one eight block) to \u2588 (full block)
//   https://en.wikipedia.org/wiki/Block_Elements
// LOWER ONE EIGHTH BLOCK … FULL BLOCK
static const char *const BLOCK_CHARS8[8] = {"\u2581", "\u2582", "\u2583", "\u2584", "\u2585", "\u2586", "\u2587", "\u2588"};
// In some terminals / with some fonts not all UTF8 block characters are supported (e.g. Putty). Use only half and full block for those
static const char *const BLOCK_CHARS2[2] = {"\u2584", "\u2588"};

#else

// Fall back to explicit UTF-8 encodings of the same characters
#define VERTICAL_LINE "\xE2\x94\x82"

static const char *const BLOCK_CHARS8[8] = {
    "\xE2\x96\x81", "\xE2\x96\x82", "\xE2\x96\x83", "\xE2\x96\x84",
    "\xE2\x96\x85", "\xE2\x96\x86", "\xE2\x96\x87", "\xE2\x96\x88" };

static const char *const BLOCK_CHARS2[2] = {"\xE2\x96\x84", "\xE2\x96\x88"};

#endif

// in bam_plcmd.c
int read_file_list(const char *file_list, int *n, char **argv[]);

static int usage() {
    fprintf(stdout, "Usage: samtools coverage [options] in1.bam [in2.bam [...]]\n\n"
            "Input options:\n"
            "  -b, --bam-list FILE     list of input BAM filenames, one per line\n"
            "  -l, --min-read-len INT  ignore reads shorter than INT bp [0]\n"
            "  -q, --min-MQ INT        base quality threshold [0]\n"
            "  -Q, --min-BQ INT        mapping quality threshold [0]\n"
            "  --rf <int|str>          required flags: skip reads with mask bits unset []\n"
            "  --ff <int|str>          filter flags: skip reads with mask bits set \n"
            "                                      [UNMAP,SECONDARY,QCFAIL,DUP]\n"
            "Output options:\n"
            "  -m, --histogram         show histogram instead of tabular output\n"
            "  -A, --ascii             show only ASCII characters in histogram\n"
            "  -o, --output FILE       write output to FILE [stdout]\n"
            "  -H, --no-header         don't print a header in tabular mode\n"
            "  -w, --n-bins INT        number of bins in histogram [terminal width - 40]\n"
            "  -r, --region REG        show specified region. Format: chr:start-end. \n"
            "  -h, --help              help (this page)\n");

    fprintf(stdout, "\nGeneric options:\n");
    sam_global_opt_help(stdout, "-.--.--.");

    fprintf(stdout,
            "\nSee manpage for additional details.\n"
            "  rname       Reference name / chromosome\n"
            "  startpos    Start position\n"
            "  endpos      End position (or sequence length)\n"
            "  numreads    Number reads aligned to the region (after filtering)\n"
            "  covbases    Number of covered bases with depth >= 1\n"
            "  coverage    Proportion of covered bases [0..1]\n"
            "  meandepth   Mean depth of coverage\n"
            "  meanbaseq   Mean baseQ in covered region\n"
            "  meanmapq    Mean mapQ of selected reads\n"
           );

    return EXIT_SUCCESS;
}

static char* center_text(char *text, char *buf, int width) {
    int len = strlen(text);
    assert(len <= width);
    int padding = (width - len) / 2;
    int padding_ex = (width - len) % 2;
    if (padding >= 1)
        sprintf(buf, " %*s%*s", len+padding, text, padding-1+padding_ex, " ");
    else
        sprintf(buf, "%s", text);

    return buf;
}

static char* readable_bps(double base_pairs, char *buf) {
    const char* units[] = {"", "K", "M", "G", "T"};
    int i = 0;
    while (base_pairs >= 1000 && i < (sizeof(units)/sizeof(units[0]) - 1)) {
        base_pairs /= 1000;
        i++;
    }
    sprintf(buf, "%.*f%s", i, base_pairs, units[i]);
    return buf;
}

static void set_read_counts(bam_aux_t **data, stats_aux_t *stats, int n_bam_files) {
    int i;
    stats->n_reads = 0;
    stats->n_selected_reads = 0;
    stats->summed_mapQ = 0;
    for (i = 0; i < n_bam_files && data[i]; ++i) {
        stats->n_reads += data[i]->n_reads;
        stats->n_selected_reads += data[i]->n_selected_reads;
        stats->summed_mapQ += data[i]->summed_mapQ;
        data[i]->n_reads = 0;
        data[i]->n_selected_reads = 0;
        data[i]->summed_mapQ = 0;
    }
}

// read one alignment from one BAM file
static int read_bam(void *data, bam1_t *b) {
    bam_aux_t *aux = (bam_aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1) {
        if((ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b)) < 0) break;
        ++aux->n_reads;

        if ( aux->fail_flags && (b->core.flag & aux->fail_flags) ) continue;
        if ( aux->required_flags && !(b->core.flag & aux->required_flags) ) continue;
        if ( b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        ++aux->n_selected_reads;
        aux->summed_mapQ += b->core.qual;
        break;
    }
    return ret;
}

void print_tabular_line(FILE *file_out, const sam_hdr_t *h, const stats_aux_t *stats) {
    fputs(sam_hdr_tid2name(h, stats->tid), file_out);
    double region_len = (double) stats->end - stats->beg;
    fprintf(file_out, "\t%"PRId64"\t%"PRId64"\t%u\t%llu\t%g\t%g\t%.3g\t%.3g\n",
            stats->beg+1,
            stats->end,
            stats->n_selected_reads,
            stats->n_covered_bases,
            100.0 * stats->n_covered_bases / region_len,
            stats->summed_coverage / region_len,
            stats->summed_coverage > 0? stats->summed_baseQ/(double) stats->summed_coverage : 0,
            stats->n_selected_reads > 0? stats->summed_mapQ/(double) stats->n_selected_reads : 0
           );
}

void print_hist(FILE *file_out, const sam_hdr_t *h, const stats_aux_t *stats, const uint32_t *hist,
        const int hist_size, const bool full_utf) {
    int i, col;
    bool show_percentiles = false;
    const int n_rows = 10;
    const char * const * BLOCK_CHARS = full_utf? BLOCK_CHARS8 : BLOCK_CHARS2;
    const int blockchar_len = full_utf? 8 : 2;
    /*
       if (stats->beg == 0) {
       stats->end = h->target_len[stats->tid];
       }
       */
    double region_len = stats->end - stats->beg;

    // Calculate histogram that contains percent covered
    double hist_data[hist_size];
    double max_val = 0.0;
    for (i = 0; i < hist_size; ++i) {
        hist_data[i] = 100 * hist[i] / (double) stats->bin_width;
        if (hist_data[i] > max_val) max_val = hist_data[i];
    }

    char buf[30];
    fprintf(file_out, "%s (%sbp)\n", sam_hdr_tid2name(h, stats->tid), readable_bps(sam_hdr_tid2len(h, stats->tid), buf));

    double row_bin_size = max_val / (double) n_rows;
    for (i = n_rows-1; i >= 0; --i) {
        double current_bin = row_bin_size * i;
        if (show_percentiles) {
            fprintf(file_out, ">%3i%% ", i*10);
        } else {
            fprintf(file_out, ">%7.2f%% ", current_bin);
        }
        fprintf(file_out, VERTICAL_LINE);
        for (col = 0; col < hist_size; ++col) {
            // get the difference in eights, or halfs when full UTF8 is not supported
            int cur_val_diff = round(blockchar_len * (hist_data[col] - current_bin) / row_bin_size) - 1;
            if (cur_val_diff < 0) {
                fputc(' ', file_out);
            } else {
                if (cur_val_diff >= blockchar_len)
                    cur_val_diff = blockchar_len - 1;

                fprintf(file_out, "%s", BLOCK_CHARS[cur_val_diff]);
            }
        }
        fprintf(file_out, VERTICAL_LINE);
        fputc(' ', file_out);
        switch (i) {
            case 9: fprintf(file_out, "Number of reads: %i", stats->n_selected_reads); break;
            case 8: if (stats->n_reads - stats->n_selected_reads > 0) fprintf(file_out, "    (%i filtered)", stats->n_reads - stats->n_selected_reads); break;
            case 7: fprintf(file_out, "Covered bases:   %sbp", readable_bps(stats->n_covered_bases, buf)); break;
            case 6: fprintf(file_out, "Percent covered: %.4g%%",
                            100.0 * stats->n_covered_bases / region_len); break;
            case 5: fprintf(file_out, "Mean coverage:   %.3gx",
                            stats->summed_coverage / region_len); break;
            case 4: fprintf(file_out, "Mean baseQ:      %.3g",
                            stats->summed_baseQ/(double) stats->summed_coverage); break;
            case 3: fprintf(file_out, "Mean mapQ:       %.3g",
                            stats->summed_mapQ/(double) stats->n_selected_reads); break;
            case 1: fprintf(file_out, "Histo bin width: %sbp",
                            readable_bps(stats->bin_width, buf)); break;
            case 0: fprintf(file_out, "Histo max bin:   %.5g%%", max_val); break;
        };
        fputc('\n', file_out);
    }

    // print x axis. Could be made pretty for widths that are not divisible
    // by 10 by variable spacing of the labels, instead of placing a label every 10 characters
    char buf2[50];
    fprintf(file_out, "     %s", center_text(readable_bps(stats->beg + 1, buf), buf2, 10));
    int rest;
    for (rest = 10; rest < 10*(hist_size/10); rest += 10) {
        fprintf(file_out, "%s", center_text(readable_bps(stats->beg + stats->bin_width*rest, buf), buf2, 10));
    }
    int last_padding = hist_size%10;
    fprintf(file_out, "%*s%s", last_padding, " ", center_text(readable_bps(stats->end, buf), buf2, 10));
    fprintf(file_out, "\n");
}

int main_coverage(int argc, char *argv[]) {
    int status = EXIT_SUCCESS;

    int ret, tid, pos, i, j;

    int max_depth = 0;
    int opt_min_baseQ = 0;
    int opt_min_mapQ = 0;
    int opt_min_len = 0;
    int opt_n_bins = 50;
    bool opt_full_width = true;
    char *opt_output_file = NULL;
    bam_aux_t **data = NULL;
    bam_mplp_t mplp = NULL;
    const bam_pileup1_t **plp = NULL;
    uint32_t *hist = NULL;
    stats_aux_t *stats = NULL;
    char *opt_reg = 0; // specified region
    char *opt_file_list = NULL;
    int n_bam_files = 0;
    char **fn = NULL;
    int fail_flags = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP); // Default fail flags
    int required_flags = 0;

    int *n_plp = NULL;
    sam_hdr_t *h = NULL; // BAM header of the 1st input

    bool opt_print_header = true;
    bool opt_print_tabular = true;
    bool opt_print_histogram = false;
    bool *covered_tids = NULL;
    bool opt_full_utf = true;

    FILE *file_out = stdout;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-'),
        {"rf", required_argument, NULL, 1}, // require flag
        {"ff", required_argument, NULL, 2}, // filter flag
        {"incl-flags", required_argument, NULL, 1}, // require flag
        {"excl-flags", required_argument, NULL, 2}, // filter flag
        {"bam-list", required_argument, NULL, 'b'},
        {"min-read-len", required_argument, NULL, 'L'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-mq", required_argument, NULL, 'q'},
        {"min-BQ", required_argument, NULL, 'Q'},
        {"min-bq", required_argument, NULL, 'Q'},
        {"histogram", no_argument, NULL, 'm'},
        {"ascii", no_argument, NULL, 'A'},
        {"output", required_argument, NULL, 'o'},
        {"no-header", no_argument, NULL, 'H'},
        {"n-bins", required_argument, NULL, 'w'},
        {"region", required_argument, NULL, 'r'},
        {"help", no_argument, NULL, 'h'},
        { NULL, 0, NULL, 0 }
    };

    // parse the command line
    int c;
    opterr = 0;
    while ((c = getopt_long(argc, argv, "Ao:L:q:Q:hHw:r:b:m", lopts, NULL)) != -1) {
        switch (c) {
            case 1:
                if ((required_flags = bam_str2flag(optarg)) < 0) {
                    fprintf(stderr,"Could not parse --rf %s\n", optarg); return EXIT_FAILURE;
                }; break;
            case 2:
                if ((fail_flags = bam_str2flag(optarg)) < 0) {
                    fprintf(stderr,"Could not parse --ff %s\n", optarg); return EXIT_FAILURE;
                }; break;
            case 'o': opt_output_file = optarg; opt_full_width = false; break;
            case 'L': opt_min_len = atoi(optarg); break;
            case 'q': opt_min_baseQ = atoi(optarg); break;
            case 'Q': opt_min_mapQ = atoi(optarg); break;
            case 'w': opt_n_bins = atoi(optarg); opt_full_width = false;
                      opt_print_histogram = true; opt_print_tabular = false;
                      break;
            case 'r': opt_reg = optarg; break;   // parsing a region requires a BAM header (strdup unnecessary)
            case 'b': opt_file_list = optarg; break;
            case 'm': opt_print_histogram = true; opt_print_tabular = false; break;
            case 'A': opt_full_utf = false;
                      opt_print_histogram = true; opt_print_tabular = false;
                      break;
            case 'H': opt_print_header = false; break;
            case 'h': return usage();
            default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                          /* else fall-through */
            case '?':
                if (optopt != '?') {  // '-?' appeared on command line
                    if (optopt) { // Bad short option
                        print_error("coverage", "invalid option -- '%c'", optopt);
                    } else { // Bad long option
                        // Do our best.  There is no good solution to finding
                        // out what the bad option was.
                        // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
                        if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
                            print_error("coverage", "unrecognised option '%s'",
                                        argv[optind - 1]);
                        }
                    }
                }
                return usage();
        }
    }
    if (optind == argc && !opt_file_list)
        return usage();

    // output file provided by user
    if (opt_output_file != NULL && strcmp(opt_output_file,"-")!=0) {
        file_out = fopen( opt_output_file, "w" );
        if (file_out == NULL) {
            print_error_errno("coverage", "Cannot open \"%s\" for writing.", opt_output_file);
            return EXIT_FAILURE;
        }
    }

    if (opt_n_bins <= 0 || opt_full_width) {
        // get number of columns of terminal
        const char* env_columns = getenv("COLUMNS");
        int columns = 0;
        if (env_columns == NULL) {
#ifdef _WIN32
            CONSOLE_SCREEN_BUFFER_INFO csbi;
            if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi)) {
                columns = csbi.srWindow.Right - csbi.srWindow.Left + 1;
            }
#else
            struct winsize w;
            if (ioctl(2, TIOCGWINSZ, &w) == 0)
                columns = w.ws_col;
#endif
        } else {
            columns = atoi(env_columns); // atoi(NULL) returns 0
        }

        if (columns > 60) {
            opt_n_bins = columns - 40;
        } else {
            opt_n_bins = 40;
        }
    }

    // setvbuf(file_out, NULL, _IONBF, 0); //turn off buffering

    // Open all BAM files
    if (opt_file_list) {
        // Read file names from opt_file_list into argv, and record the number of files in n_bam_files
        if (read_file_list(opt_file_list, &n_bam_files, &fn)) {
            print_error_errno("coverage", "Cannot open file list \"%s\".", opt_file_list);
            return EXIT_FAILURE;
        }
        argv = fn;
        optind = 0;
    } else {
        n_bam_files = argc - optind; // the number of BAMs on the command line
    }

    data = (bam_aux_t **)calloc(n_bam_files, sizeof(bam_aux_t*)); // data[i] for the i-th BAM file
    if (!data) {
        print_error("coverage", "Failed to allocate memory");
        status = EXIT_FAILURE;
        goto coverage_end;
    }

    for (i = 0; i < n_bam_files; ++i) {
        int rf;
        data[i] = (bam_aux_t *) calloc(1, sizeof(bam_aux_t));
        if (!data[i]) {
            print_error("coverage", "Failed to allocate memory");
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        data[i]->fp = sam_open_format(argv[optind+i], "r", &ga.in); // open BAM

        if (data[i]->fp == NULL) {
            print_error_errno("coverage", "Could not open \"%s\"", argv[optind+i]);
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
        if (opt_min_baseQ) rf |= SAM_QUAL;

        // Set CRAM options on file handle - returns 0 on success
        if (hts_set_opt(data[i]->fp, CRAM_OPT_REQUIRED_FIELDS, rf)) {
            print_error_errno("coverage", "Failed to set CRAM_OPT_REQUIRED_FIELDS value");
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            print_error_errno("coverage", "Failed to set CRAM_OPT_DECODE_MD value");
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        data[i]->min_mapQ = opt_min_mapQ;            // set the mapQ filter
        data[i]->min_len  = opt_min_len;             // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        data[i]->fail_flags = fail_flags;
        data[i]->required_flags = required_flags;
        if (data[i]->hdr == NULL) {
            print_error_errno("coverage", "Could not read header for \"%s\"", argv[optind+i]);
            status = EXIT_FAILURE;
            goto coverage_end;
        }

        // Lookup region if specified
        if (opt_reg) { // if a region is specified
            hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
            if (idx == NULL) {
                print_error_errno("coverage", "Failed to load index for \"%s\"", argv[optind+i]);
                status = EXIT_FAILURE;
                goto coverage_end;
            }
            data[i]->iter = sam_itr_querys(idx, data[i]->hdr, opt_reg); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
            if (data[i]->iter == NULL) {
                print_error_errno("coverage", "Failed to parse region \"%s\"", opt_reg);
                status = EXIT_FAILURE;
                goto coverage_end;
            }
        }
    }

    if (opt_print_tabular && opt_print_header)
        fputs("#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n", file_out);

    h = data[0]->hdr; // easy access to the header of the 1st BAM
    int n_targets = sam_hdr_nref(h);
    covered_tids = calloc(n_targets, sizeof(bool));
    stats = calloc(1, sizeof(stats_aux_t));
    if (!covered_tids || !stats) {
        print_error("coverage", "Failed to allocate memory");
        status = EXIT_FAILURE;
        goto coverage_end;
    }

    int64_t n_bins = opt_n_bins;
    if (opt_reg) {
        stats->tid = data[0]->iter->tid;
        stats->beg = data[0]->iter->beg; // and to the parsed region coordinates
        stats->end = data[0]->iter->end;
        if (stats->end == HTS_POS_MAX) {
            stats->end = sam_hdr_tid2len(h, stats->tid);
        }
        if (opt_n_bins > stats->end - stats->beg) {
            n_bins = stats->end - stats->beg;
        }
        stats->bin_width = (stats->end-stats->beg) / n_bins;
    } else {
        stats->tid = -1;
    }

    int64_t current_bin = 0;

    // the core multi-pileup loop
    mplp = bam_mplp_init(n_bam_files, read_bam, (void**)data); // initialization
    if (max_depth > 0)
        bam_mplp_set_maxcnt(mplp, max_depth);  // set maximum coverage depth
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp, INT_MAX);


    // Extra info for histogram and coverage counting
    hist = (uint32_t*) calloc(opt_n_bins, sizeof(uint32_t));
    n_plp = (int*) calloc(n_bam_files, sizeof(int*)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = (const bam_pileup1_t**) calloc(n_bam_files, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    if (!hist || !n_plp || !plp) {
        print_error("coverage", "Failed to allocate memory");
        status = EXIT_FAILURE;
        goto coverage_end;
    }
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position

        if (tid != stats->tid) { // Next target sequence
            if (stats->tid >= 0) { // It's not the first sequence, print results
                set_read_counts(data, stats, n_bam_files);
                if (opt_print_histogram) {
                    print_hist(file_out, h, stats, hist, n_bins, opt_full_utf);
                    fputc('\n', file_out);
                } else if (opt_print_tabular) {
                    print_tabular_line(file_out, h, stats);
                }

                // reset data
                memset(stats, 0, sizeof(stats_aux_t));
                if (opt_print_histogram)
                    memset(hist, 0, n_bins*sizeof(uint32_t));
            }

            stats->tid = tid;
            covered_tids[tid] = true;
            if (!opt_reg)
                stats->end = sam_hdr_tid2len(h, tid);

            if (opt_print_histogram) {
                n_bins = opt_n_bins > stats->end-stats->beg? stats->end-stats->beg : opt_n_bins;
                stats->bin_width = (stats->end-stats->beg) / n_bins;
            }
        }
        if (pos < stats->beg || pos >= stats->end) continue; // out of range; skip
        if (tid >= n_targets) continue;     // diff number of @SQ lines per file?

        if (opt_print_histogram) {
            current_bin = (pos - stats->beg) / stats->bin_width;
        }

        bool count_base = false;
        for (i = 0; i < n_bam_files; ++i) { // base level filters have to go here
            int depth_at_pos = n_plp[i];
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modify plp[][] unless you really know

                if (p->is_del || p->is_refskip) --depth_at_pos; // having dels or refskips at tid:pos
                else if (p->qpos < p->b->core.l_qseq &&
                        bam_get_qual(p->b)[p->qpos] < opt_min_baseQ) --depth_at_pos; // low base quality
                else
                    stats->summed_baseQ += bam_get_qual(p->b)[p->qpos];
            }
            if (depth_at_pos > 0) {
                count_base = true;
                stats->summed_coverage += depth_at_pos;
            }
            // hist[current_bin] += depth_at_pos;  // Add counts to the histogram here to have one based on coverage
            //fprintf(file_out, "\t%d", n_plp[i] - m); // this the depth to output
        }
        if (count_base) {
            ++(stats->n_covered_bases);
            if (opt_print_histogram && current_bin < n_bins)
                ++(hist[current_bin]); // Histogram based on breadth of coverage
        }
    }

    if (stats->tid != -1) {
        set_read_counts(data, stats, n_bam_files);
        if (opt_print_histogram) {
            print_hist(file_out, h, stats, hist, n_bins, opt_full_utf);
        } else if (opt_print_tabular) {
            print_tabular_line(file_out, h, stats);
        }
    }


    if (!opt_reg && opt_print_tabular) {
        memset(stats, 0, sizeof(stats_aux_t));
        for (i = 0; i < n_targets; ++i) {
            if (!covered_tids[i]) {
                stats->tid = i;
                stats->end = sam_hdr_tid2len(h, i);
                print_tabular_line(file_out, h, stats);
            }
        }
    }

    if (ret < 0) status = EXIT_FAILURE;

coverage_end:
    if (n_plp) free(n_plp);
    if (plp) free(plp);
    bam_mplp_destroy(mplp);

    if (covered_tids) free(covered_tids);
    if (hist) free(hist);
    if (stats) free(stats);


    // Close files and free data structures
    if (!(file_out == stdout || fclose(file_out) == 0)) {
        if (status == EXIT_SUCCESS) {
            print_error_errno("coverage", "error on closing \"%s\"",
                    (opt_output_file && strcmp(opt_output_file, "-") != 0?
                     opt_output_file : "stdout"));
            status = EXIT_FAILURE;
        }
    }

    if (data) {
        for (i = 0; i < n_bam_files && data[i]; ++i) {
            sam_hdr_destroy(data[i]->hdr);
            if (data[i]->fp) sam_close(data[i]->fp);
            hts_itr_destroy(data[i]->iter);
            free(data[i]);
        }
        free(data);
    }

    if (opt_file_list && fn) {
        for (i = 0; i < n_bam_files; ++i)
            free(fn[i]);
        free(fn);
    }
    sam_global_args_free(&ga);

    return status;
}

#ifdef _MAIN_BAMCOV
int main(int argc, char *argv[]) {
    return main_coverage(argc, argv);
}
#endif
