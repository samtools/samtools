/* coverage.c -- samtools coverage subcommand
 
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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>  // variadic functions
#include <limits.h>  // INT_MAX
#include <math.h>    // round
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <sys/ioctl.h>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "samtools.h"
#include "sam_opts.h"

const char *VERSION = "0.1";

typedef struct {  // auxiliary data structure to hold a BAM file
    samFile *fp;     // file handle
    bam_hdr_t *hdr;  // file header
    hts_itr_t *iter; // iterator to a region - NULL for us by default
    int min_mapQ;    // mapQ filter
    int min_len;     // length filter
    uint32_t n_reads;  // records the number of reads seen in file
    uint32_t n_selected_reads; // records the number of reads passing filter
    uint32_t summed_mapQ; // summed mapQ of all reads passing filter
} bam_aux_t;

typedef struct {  // auxiliary data structure to hold stats on coverage
    int32_t tid;    // chromosome ID, defined by header
    int beg;
    int end;
    int bin_width;
    uint64_t n_covered_bases;
    uint64_t summed_coverage;
    uint64_t summed_baseQ;
    uint64_t summed_mapQ;
    uint32_t n_reads;
    uint32_t n_selected_reads;
} stats_aux_t;

int read_file_list(const char *file_list, int *n, char **argv[]);

static int usage(int status) {
    fprintf(stderr, "\nUsage: samtools coverage [options] in1.bam [in2.bam [...]]\n"
                    "Options:\n"
                    "  -o FILE           where to write output to [stdout]\n"
                    "  -l <int>          read length threshold - ignore reads shorter than <int> [0]\n"
                    "  -q <int>          base quality threshold [0]\n"
                    "  -Q <int>          mapping quality threshold [0]\n"
                    "  -H                print a file header\n"
                    "  -m                show histogram instead of tabular output\n"
                    "  -w <int>          number of bins in histogram. Set to 0 to use terminal width [50]\n"
                    "  -r <chr:from-to>  Show region on chromosome chr. \n"
                    "  -h                help (this page)\n"
                    "  -v                version of this command\n"
                    "\nGlobal options:\n");

    sam_global_opt_help(stderr, "-.--.-");

    fprintf(stderr, 
            "\nThe tabular output is a simple tab-seperate table with 7 columns:\n"
                    "  1. Reference name / chromosome\n"
                    "  2. Start position\n"
                    "  3. End position (or sequence length)\n"
                    "  4. Number reads aligned to the region\n"
                    "  5. Number of covered bases with depth >= 1\n"
                    "  6. Percent of covered bases\n"
                    "  7. Average coverage depth\n"
                    "  8. Average baseQ in covered region\n"
                    "  9. Average mapQ of selected reads\n"
                    );

    return status;
}

static int version(int status) {
    fprintf(stderr,
    "┌┐ ┌─┐┌┬┐┌─┐┌─┐┬  ┬\n"
    "├┴┐├─┤││││  │ │└┐┌┘\n"
    "└─┘┴ ┴┴ ┴└─┘└─┘ └┘   v%s\n", VERSION);
    return status;
}

char* centerText(char *text, char *buf, int width) {
    int len = strlen(text);
    int padding = (width - len) / 2;
    int padding_ex = (width - len) % 2;
    if (padding >= 1) 
        sprintf(buf, " %*s%*s", len+padding, text, padding-1+padding_ex, " ");
    return(buf);
}

char* readable_bps(double base_pairs, char *buf) {
    const char* units[] = {"", "K", "M", "G", "T"};
    int i = 0;
    while (base_pairs >= 1000) {
        base_pairs /= 1000;
        i++;
    }
    sprintf(buf, "%.*f%s", i, base_pairs, units[i]);
    return buf;
}

void set_read_counts(bam_aux_t **data, stats_aux_t *stats, int n_bam_files) {
    stats->n_reads = 0;
    stats->n_selected_reads = 0;
    stats->summed_mapQ = 0;
    for (int i = 0; i < n_bam_files && data[i]; ++i) {
        stats->n_reads += data[i]->n_reads;
        stats->n_selected_reads += data[i]->n_selected_reads;
        stats->summed_mapQ += data[i]->summed_mapQ;
        data[i]->n_reads = 0;
        data[i]->n_selected_reads = 0;
        data[i]->summed_mapQ = 0;
    }
}

// UTF8 specifies block characters in eights going from \u2581 (lower one eight block) to \u2588 (full block)
//   https://en.wikipedia.org/wiki/Block_Elements
const char* BLOCK_CHARS[8] = {"\u2581", "\u2582", "\u2583", "\u2584", "\u2585", "\u2586", "\u2587", "\u2588"};

// read one alignment from one BAM file
static int read_bam(void *data, bam1_t *b) {
    bam_aux_t *aux = (bam_aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1) {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        ++aux->n_reads;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        ++aux->n_selected_reads;
        aux->summed_mapQ += b->core.qual;
        break;
    }
    return ret;
}

void print_tabular_line(FILE *file_out, const bam_hdr_t *h, const stats_aux_t *stats) {
    fputs(h->target_name[stats->tid], file_out);
    double region_len = (double) stats->end - stats->beg;
    fprintf(file_out, "\t%d\t%d\t%u\t%llu\t%g\t%g\t%.3g\t%.3g\n",
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

void print_hist(FILE *file_out, const bam_hdr_t *h, const stats_aux_t *stats, const uint32_t *hist, 
                const int hist_size) {

    bool show_percentiles = false;
    const int n_rows = 10;
    /*
    if (stats->beg == 0) {
        stats->end = h->target_len[stats->tid];
    }
    */
    double region_len = stats->end - stats->beg;

    // Calculate histogram that contains percent covered
    double hist_d[hist_size];
    for (int i = 0; i < hist_size; ++i) {
        hist_d[i] = 100 * hist[i] / (double) stats->bin_width;
        //std::cerr << col_bin_size * i << "," << hist_d[i] << "\n";
        //std::cerr << "{x:" << col_bin_size * i << ",y:" << hist_d[i] << "},";
    }
    //std::cerr << std::endl;
    //return;

    double max_val = 0.0;
    for (int i = 0; i < hist_size; ++i) {
        if (hist_d[i] > max_val) max_val = hist_d[i];
    }

    char buf[100];

    fprintf(file_out, "%s (%sbp)\n", h->target_name[stats->tid], readable_bps(h->target_len[stats->tid], buf));
    
    double row_bin_size = max_val / (double) n_rows;

    for (int row = n_rows-1; row >= 0; --row) {
        double current_bin = row_bin_size * row;
        if (show_percentiles) {
            fprintf(file_out, ">%3i%% \u2502", row*10);
        } else {
            //if (row == 0) 
            //    fprintf(file_out, ">%7.2f%% \u23A3", current_bin);
            //else
            fprintf(file_out, ">%7.2f%% \u2502", current_bin);
        }
        for (int col = 0; col < hist_size; ++col) {
            // get the difference in eights
            int cur_val_diff = round(8 * (hist_d[col] - current_bin) / row_bin_size);
            if (cur_val_diff == 0) {
                fputc('_', file_out);
            } else if (cur_val_diff > 0) {
                if (cur_val_diff > 8)
                    cur_val_diff = 8;
                fprintf(file_out, "%s", BLOCK_CHARS[cur_val_diff-1]);
            } else {
                fputc(' ', file_out);
            }
        }
        fprintf(file_out, "\u2502 ");
        //fprintf(file_out, "  ");
        switch (row) {
            case 9: fprintf(file_out, "Number of reads:  %i", stats->n_selected_reads); break;
            case 8: if (stats->n_reads - stats->n_selected_reads > 0) fprintf(file_out, "    (%i filtered)", stats->n_reads - stats->n_selected_reads); break;
            case 7: fprintf(file_out, "Covered bases:    %sbp", readable_bps(stats->n_covered_bases, buf)); break;
            case 6: fprintf(file_out, "Percent covered:  %.4g%%", 
                            100.0 * stats->n_covered_bases / region_len); break;
            case 5: fprintf(file_out, "Average coverage: %.3gx",
                            stats->summed_coverage / region_len); break;
            case 4: fprintf(file_out, "Average baseQ:    %.3g",
                            stats->summed_baseQ/(double) stats->summed_coverage); break;
            case 3: fprintf(file_out, "Average mapQ:     %.3g",
                            stats->summed_mapQ/(double) stats->n_selected_reads); break;
            case 1: fprintf(file_out, "Histo bin width:  %sbp",
                            readable_bps(stats->bin_width, buf)); break;
            case 0: fprintf(file_out, "Histo max bin:    %.5g%%", max_val); break;
        };
        fputc('\n', file_out);
    }

    // print x axis. Could be made pretty for widths that are not divisible 
    // by 10 by variable spacing of the labels, instead of placing a label every 10 characters
    char buf2[11];
    fprintf(file_out, "     %s", centerText(readable_bps(stats->beg + 1, buf), buf2, 10));
    for (int rest = 10; rest < 10*(hist_size/10); rest += 10) {
        fprintf(file_out, "%s", centerText(readable_bps(stats->beg + stats->bin_width*rest, buf), buf2, 10));
    }
    int last_padding = hist_size%10;
    fprintf(file_out, "%*s%s", last_padding, " ", centerText(readable_bps(stats->end, buf), buf2, 10));
    fprintf(file_out, "\n");
}

int main_coverage(int argc, char *argv[]) {
    int status = EXIT_SUCCESS;

    int ret, tid, pos, i;

    int max_depth = 0;
    int opt_min_baseQ = 0;
    int opt_min_mapQ = 0;
    int opt_min_len = 0;
    int opt_n_bins = 50;
    char *opt_output_file = NULL;
    bam_aux_t **data;
    bam_mplp_t mplp;
    const bam_pileup1_t **plp;
    char *opt_reg = 0; // specified region
    char *opt_file_list = NULL;
    int n_bam_files = 0;

    int *n_plp;
    bam_hdr_t *h = NULL; // BAM header of the 1st input

    bool opt_print_header = true;
    bool opt_print_tabular = true;
    bool opt_print_histogram = false;
    bool opt_print_html = false;
    bool *covered_tids;

    int current_tid = -1;
    FILE *file_out = stdout;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-'),
        { NULL, 0, NULL, 0 }
    };

    // parse the command line
    int c;
    while ((c = getopt_long(argc, argv, "o:l:q:Q:hH:w:vr:Df:m", lopts, NULL)) != -1) {
            switch (c) {
                case 'o': opt_output_file = optarg; break;
                case 'l': opt_min_len = atoi(optarg); break;
                case 'q': opt_min_baseQ = atoi(optarg); break;
                case 'Q': opt_min_mapQ = atoi(optarg); break;
                case 'w': opt_n_bins = atoi(optarg); break;
                case 'r': opt_reg = strdup(optarg); break;   // parsing a region requires a BAM header
                case 'f': opt_file_list = optarg; break;
                case 'm': opt_print_histogram = true; opt_print_tabular = false; break;
                case 'H': opt_print_header = true; break;
                case 'D': opt_print_html = true; opt_print_tabular = false; break;
                case 'h': return usage(EXIT_SUCCESS);
                case 'v': return version(EXIT_SUCCESS);
                default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
                case '?': return usage(EXIT_FAILURE);
            }
    }

    bool compute_histogram = opt_print_histogram || opt_print_html;

    if (optind == argc && !opt_file_list) 
        return usage(EXIT_SUCCESS);

    // output file provided by user
    if (opt_output_file != NULL && strcmp(opt_output_file,"-")!=0) {
        file_out = fopen( opt_output_file, "w" );
        if (file_out == NULL) {
            print_error_errno("coverage", "Cannot open \"%s\" for writing.", opt_output_file);
            return EXIT_FAILURE;
        }
    }

    if (opt_n_bins <= 0) {
        // get number of lines and columns of terminal
        struct winsize w; // lines are in field ws_row, columns in field ws_col
        if (ioctl(2, TIOCGWINSZ, &w) == 0 && w.ws_col > 60) {
            opt_n_bins = w.ws_col - 40;
        }
    }

    //setvbuf(file_out, NULL, _IONBF, 0); //turn off buffering

    // Open all BAM files
    
    if (opt_file_list) {
        // Read file names from opt_file_list into argv, and record the number of files in n_bam_files
        if ( read_file_list(opt_file_list, &n_bam_files, &argv) ) {
            print_error_errno("coverage", "Cannot open file list \"%s\".", opt_file_list);
            return EXIT_FAILURE;
        }
        optind = 0;
    }
    else
        n_bam_files = argc - optind; // the number of BAMs on the command line

    data = (bam_aux_t **)calloc(n_bam_files, sizeof(bam_aux_t*)); // data[i] for the i-th BAM file

    for (int i = 0; i < n_bam_files; ++i) {
        int rf;
        data[i] = (bam_aux_t *) calloc(1, sizeof(bam_aux_t));
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
        if (data[i]->hdr == NULL) {
            print_error_errno("coverage", "Couldn't read header for \"%s\"", argv[optind+i]);
            status = EXIT_FAILURE;
            goto coverage_end;
        }

        // Lookup region if specified
        if (opt_reg) { // if a region is specified
            hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
            if (idx == NULL) {
                print_error_errno("coverage", "can't load index for \"%s\"", argv[optind+i]);
                status = EXIT_FAILURE;
                goto coverage_end;
            }
            data[i]->iter = sam_itr_querys(idx, data[i]->hdr, opt_reg); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
            if (data[i]->iter == NULL) {
                print_error_errno("coverage", "can't parse region \"%s\"", opt_reg);
                status = EXIT_FAILURE;
                goto coverage_end;
            }
        }
    }

    if (opt_print_tabular && opt_print_header)
        fputs("#CHROM\tSTART\tEND\tN_READS\tN_COVERED_BASES\tPERCENT_COVERED\tAVG_COV\tAVG_BASEQ\tAVG_MAPQ\n", file_out);


    h = data[0]->hdr; // easy access to the header of the 1st BAM
    covered_tids = calloc(h->n_targets, sizeof(bool));
    stats_aux_t *stats = calloc(1, sizeof(stats_aux_t));

    bool reg_sets_end = false;
    if (opt_reg) {
        stats->beg = data[0]->iter->beg; // and to the parsed region coordinates
        stats->end = data[0]->iter->end;
        stats->tid = data[0]->iter->tid;
        stats->bin_width = (stats->end-stats->beg) / opt_n_bins;
        reg_sets_end = stats->end != INT_MAX;
    } else {
        stats->tid = -1;
    }

    int current_bin;

    // the core multi-pileup loop
    mplp = bam_mplp_init(n_bam_files, read_bam, (void**)data); // initialization
    if (max_depth > 0)
        bam_mplp_set_maxcnt(mplp,max_depth);  // set maximum coverage depth
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp,INT_MAX);


    // Extra info for histogram and coverage counting
    uint32_t *hist = (uint32_t*) calloc(opt_n_bins, sizeof(uint32_t));
    //uint32_t **hists_for_html = (uint32_t**) calloc(opt_n_bins*h->n_targets, sizeof(uint32_t));
    //stats_aux_t **stats_for_html = (stats_aux_t*) calloc(h->n_targets, sizeof(stats_aux_t));
    
    n_plp = (int*) calloc(n_bam_files, sizeof(int*)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = (const bam_pileup1_t**) calloc(n_bam_files, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position
        if (tid != stats->tid) {
            if (stats->tid >= 0) {
                set_read_counts(data, stats, n_bam_files);
                if (opt_print_histogram) {
                    print_hist(file_out, h, stats, hist, opt_n_bins);
                    fputc('\n', file_out);
                } else if (opt_print_tabular) {
                    print_tabular_line(file_out, h, stats);
                } else if (opt_print_html) {
                    //hists_for_html.push_back(hist);
                    //stats_for_html.push_back(stats);
                }

                for (i = 0; i < n_bam_files && data[i]; ++i) {
                    data[i]->n_reads = 0;
                    data[i]->n_selected_reads = 0;
                }

                memset(stats, 0, sizeof(stats_aux_t));

                if (compute_histogram)
                    memset(hist, 0, opt_n_bins*sizeof(uint32_t));
            }

            stats->tid = tid;
            covered_tids[tid] = true;
            if (!reg_sets_end)
                stats->end = h->target_len[tid];

            if (compute_histogram) 
                stats->bin_width = (stats->end-stats->beg) / opt_n_bins;
        }
        if (pos < stats->beg || pos >= stats->end) continue; // out of range; skip
        if (tid >= h->n_targets) continue;     // diff number of @SQ lines per file?

        if (compute_histogram) {
            current_bin = (pos - stats->beg) / stats->bin_width;
        }

        bool count_base = false;
        for (i = 0; i < n_bam_files; ++i) { // base level filters have to go here
            int depth_at_pos = n_plp[i];
            for (int j = 0; j < n_plp[i]; ++j) {
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
            if (compute_histogram && current_bin < opt_n_bins)
                ++(hist[current_bin]); // Histogram based on breadth of coverage
        }
    }

    if (stats->tid != -1) {
        set_read_counts(data, stats, n_bam_files);
        if (opt_print_histogram) {
            print_hist(file_out, h, stats, hist, opt_n_bins);
        } else if (opt_print_tabular) {
            print_tabular_line(file_out, h, stats);
        }
    }


    if (!opt_reg && opt_print_tabular) {
        memset(stats, 0, sizeof(stats_aux_t));
        for (int tid = 0; tid < h->n_targets; ++tid) {
            if (!covered_tids[tid]) {
                stats->tid = tid;
                stats->end = h->target_len[tid];
                print_tabular_line(file_out, h, stats);
            }
        }
    }

    if (ret < 0) status = EXIT_FAILURE;

    free(n_plp); 
    free(plp);
    bam_mplp_destroy(mplp);

    free(covered_tids);
    free(hist);
    free(stats);

coverage_end:
    // Close files and free data structures
    if (fclose(file_out) != 0) {
        if (status == EXIT_SUCCESS) {
            print_error_errno("coverage", "error on closing \"%s\"", 
                    (opt_output_file && strcmp(opt_output_file, "-") != 0?
                     opt_output_file : "stdout"));
            status = EXIT_FAILURE;
        }
    }
    for (i = 0; i < n_bam_files && data[i]; ++i) {
        bam_hdr_destroy(data[i]->hdr);
        if (data[i]->fp) sam_close(data[i]->fp);
        hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }

    free(data); 
    free(opt_reg);
    if (opt_file_list) {
        for (int i = 0; i < n_bam_files; ++i) {
            free(argv[i]);
            free(argv);
        }
    }
    sam_global_args_free(&ga);

    return status;
}

#ifdef _MAIN_BAMCOV
int main(int argc, char *argv[]) {
    return main_coverage(argc, argv);
}
#endif
