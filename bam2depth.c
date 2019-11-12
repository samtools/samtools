/*  bam2depth.c -- depth subcommand.

    Copyright (C) 2011, 2012 Broad Institute.
    Copyright (C) 2012-2016, 2018, 2019 Genome Research Ltd.

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

/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -lhts -lz
 */

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "samtools.h"
#include "bedidx.h"
#include "sam_opts.h"

#define BAM_FMAX ((BAM_FSUPPLEMENTARY << 1) - 1)

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    sam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
    uint32_t flags;  // read filtering flags
} aux_t;

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1)
    {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & aux->flags) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        break;
    }
    return ret;
}

int read_file_list(const char *file_list,int *n,char **argv[]);

static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: samtools depth [options] in1.bam [in2.bam [...]]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -a                  output all positions (including zero depth)\n");
    fprintf(stderr, "   -a -a (or -aa)      output absolutely all positions, including unused ref. sequences\n");
    fprintf(stderr, "   -b <bed>            list of positions or regions\n");
    fprintf(stderr, "   -X                  use customized index files\n");
    fprintf(stderr, "   -f <list>           list of input BAM filenames, one per line [null]\n");
    fprintf(stderr, "   -H                  print a file header\n");
    fprintf(stderr, "   -l <int>            read length threshold (ignore reads shorter than <int>) [0]\n");
    fprintf(stderr, "   -d/-m <int>         maximum coverage depth [8000]. If 0, depth is set to the maximum\n"
                    "                       integer value, effectively removing any depth limit.\n");  // the htslib's default
    fprintf(stderr, "   -o FILE             where to write output to [stdout]\n");
    fprintf(stderr, "   -q <int>            base quality threshold [0]\n");
    fprintf(stderr, "   -Q <int>            mapping quality threshold [0]\n");
    fprintf(stderr, "   -r <chr:from-to>    region\n");
    fprintf(stderr, "   -g <flags>          include reads that have any of the specified flags set [0]\n");
    fprintf(stderr, "   -G <flags>          filter out reads that have any of the specified flags set"
                    "                       [UNMAP,SECONDARY,QCFAIL,DUP]\n");

    sam_global_opt_help(stderr, "-.--.--.");

    fprintf(stderr, "\n");
    fprintf(stderr, "The output is a simple tab-separated table with three columns: reference name,\n");
    fprintf(stderr, "position, and coverage depth.  Note that positions with zero coverage may be\n");
    fprintf(stderr, "omitted by default; see the -a option.\n");
    fprintf(stderr, "\n");

    return EXIT_FAILURE;
}

int main_depth(int argc, char *argv[])
{
    int i, n, tid, reg_tid, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, has_index_file = 0;
    hts_pos_t beg, end, pos, last_pos = -1;
    int all = 0, status = EXIT_SUCCESS, nfiles, max_depth = -1;
    const bam_pileup1_t **plp;
    char *reg = 0; // specified region
    void *bed = 0; // BED data structure
    char *file_list = NULL, **fn = NULL;
    sam_hdr_t *h = NULL; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int last_tid = -1, ret;
    int print_header = 0;
    char *output_file = NULL;
    FILE *file_out = stdout;
    uint32_t flags = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
    int tflags = 0;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-'),
        { NULL, 0, NULL, 0 }
    };

    // parse the command line
    while ((n = getopt_long(argc, argv, "r:b:Xq:Q:l:f:am:d:Ho:g:G:", lopts, NULL)) >= 0) {
        switch (n) {
            case 'l': min_len = atoi(optarg); break; // minimum query length
            case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
            case 'b':
                bed = bed_read(optarg); // BED or position list file can be parsed now
                if (!bed) {
                    print_error_errno("depth", "Could not read file \"%s\"", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'X': has_index_file = 1; break;
            case 'q': baseQ = atoi(optarg); break;   // base quality threshold
            case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
            case 'f': file_list = optarg; break;
            case 'a': all++; break;
            case 'd': case 'm': max_depth = atoi(optarg); break; // maximum coverage depth
            case 'H': print_header = 1; break;
            case 'o': output_file = optarg; break;
            case 'g':
                tflags = bam_str2flag(optarg);
                if (tflags < 0 || tflags > BAM_FMAX) {
                    print_error_errno("depth", "Flag value \"%s\" is not supported", optarg);
                    return 1;
                }
                flags &= ~tflags;
                break;
            case 'G':
                tflags = bam_str2flag(optarg);
                if (tflags < 0 || tflags > BAM_FMAX) {
                    print_error_errno("depth", "Flag value \"%s\" is not supported", optarg);
                    return 1;
                }
                flags |= tflags;
                break;
            default:  if (parse_sam_global_opt(n, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
            case '?': return usage();
        }
    }
    if (optind == argc && !file_list)
        return usage();

    /* output file provided by user */
    if (output_file != NULL && strcmp(output_file,"-")!=0) {
        file_out = fopen( output_file, "w" );
        if (file_out == NULL) {
            print_error_errno("depth", "Cannot open \"%s\" for writing.", output_file);
            return EXIT_FAILURE;
        }
    }


    // initialize the auxiliary data structures
    if (file_list)
    {
        if (has_index_file) {
            print_error("depth", "The -f option cannot be combined with -X");
            return 1;
        }
        if ( read_file_list(file_list,&nfiles,&fn) ) return EXIT_FAILURE;
        n = nfiles;
        argv = fn;
        optind = 0;
    }
    else if (has_index_file) { // Calculate # of input BAM files
        if ((argc - optind) % 2 != 0) {
            fprintf(stderr, "Error: Odd number of filenames detected! Each BAM file should have an index file\n");
            return 1;
        }
        n = (argc - optind) / 2;
    } else {
        n = argc - optind;
    }
    data = calloc(n, sizeof(aux_t*)); // data[i] for the i-th input
    reg_tid = 0; beg = 0; end = HTS_POS_MAX;  // set the default region

    for (i = 0; i < n; ++i) {
        int rf;
        data[i] = calloc(1, sizeof(aux_t));
        data[i]->fp = sam_open_format(argv[optind+i], "r", &ga.in); // open BAM
        if (data[i]->fp == NULL) {
            print_error_errno("depth", "Could not open \"%s\"", argv[optind+i]);
            status = EXIT_FAILURE;
            goto depth_end;
        }
        rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
        if (baseQ) rf |= SAM_QUAL;
        if (hts_set_opt(data[i]->fp, CRAM_OPT_REQUIRED_FIELDS, rf)) {
            print_error_errno("depth", "Failed to set CRAM_OPT_REQUIRED_FIELDS value");
            status = EXIT_FAILURE;
            goto depth_end;
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            print_error_errno("depth", "Failed to set CRAM_OPT_DECODE_MD value");
            status = EXIT_FAILURE;
            goto depth_end;
        }
        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = min_len;                 // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        if (data[i]->hdr == NULL) {
            print_error_errno("depth", "Couldn't read header for \"%s\"",
                              argv[optind+i]);
            status = EXIT_FAILURE;
            goto depth_end;
        }
        if (reg) { // if a region is specified
            hts_idx_t *idx = NULL;
            // If index filename has not been specfied, look in BAM folder
            if (has_index_file) {
                idx = sam_index_load2(data[i]->fp, argv[optind+i], argv[optind+i+n]);  // load the index
            } else {
                idx = sam_index_load(data[i]->fp, argv[optind+i]);
            }
            if (idx == NULL) {
                print_error("depth", "can't load index for \"%s\"", argv[optind+i]);
                status = EXIT_FAILURE;
                goto depth_end;
            }
            data[i]->iter = sam_itr_querys(idx, data[i]->hdr, reg); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
            if (data[i]->iter == NULL) {
                print_error("depth", "can't parse region \"%s\"", reg);
                status = EXIT_FAILURE;
                goto depth_end;
            }
        }
        data[i]->flags = flags;
    }
    if (print_header) {
        fputs("#CHROM\tPOS", file_out);
        for (i = 0; i < n; ++i) {
            fputc('\t', file_out);
            fputs(argv[optind+i], file_out);
            }
        fputc('\n', file_out);
        }
    h = data[0]->hdr; // easy access to the header of the 1st BAM
    if (reg) {
        beg = data[0]->iter->beg; // and to the parsed region coordinates
        end = data[0]->iter->end;
        reg_tid = data[0]->iter->tid;
    }

    // the core multi-pileup loop
    mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
    if (0 < max_depth)
        bam_mplp_set_maxcnt(mplp,max_depth);  // set maximum coverage depth
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp,INT_MAX);
    n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    while ((ret=bam_mplp64_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        if (tid >= sam_hdr_nref(h)) continue;     // diff number of @SQ lines per file?
        if (all) {
            while (tid > last_tid) {
                if (last_tid >= 0 && !reg) {
                    // Deal with remainder or entirety of last tid.
                    while (++last_pos < sam_hdr_tid2len(h, last_tid)) {
                        // Horribly inefficient, but the bed API is an obfuscated black box.
                        if (bed && bed_overlap(bed, sam_hdr_tid2name(h, last_tid), last_pos, last_pos + 1) == 0)
                            continue;
                        fputs(sam_hdr_tid2name(h, last_tid), file_out);
                        fprintf(file_out, "\t%"PRIhts_pos, last_pos+1);
                        for (i = 0; i < n; i++)
                            fputc('\t', file_out), fputc('0', file_out);
                        fputc('\n', file_out);
                    }
                }
                last_tid++;
                last_pos = -1;
                if (all < 2)
                    break;
            }

            // Deal with missing portion of current tid
            while (++last_pos < pos) {
                if (last_pos < beg) continue; // out of range; skip
                if (bed && bed_overlap(bed, sam_hdr_tid2name(h, tid), last_pos, last_pos + 1) == 0)
                    continue;
                fputs(sam_hdr_tid2name(h, tid), file_out);
                fprintf(file_out, "\t%"PRIhts_pos, last_pos+1);
                for (i = 0; i < n; i++)
                    fputc('\t', file_out), fputc('0', file_out);
                fputc('\n', file_out);
            }

            last_tid = tid;
            last_pos = pos;
        }
        if (bed && bed_overlap(bed, sam_hdr_tid2name(h, tid), pos, pos + 1) == 0) continue;
        fputs(sam_hdr_tid2name(h, tid), file_out);
        fprintf(file_out, "\t%"PRIhts_pos, pos+1); // a customized printf() would be faster
        for (i = 0; i < n; ++i) { // base level filters have to go here
            int j, m = 0;
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
                if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
                else if (p->qpos < p->b->core.l_qseq &&
                         bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
            }
            fprintf(file_out, "\t%d", n_plp[i] - m); // this the depth to output
        }
        fputc('\n', file_out);
    }
    if (ret < 0) status = EXIT_FAILURE;
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);

    if (all) {
        // Handle terminating region
        if (last_tid < 0 && reg) {
            last_tid = reg_tid;
            last_pos = beg-1;
        }
        while (last_tid >= 0 && last_tid < sam_hdr_nref(h)) {
            while (++last_pos < sam_hdr_tid2len(h, last_tid)) {
                if (last_pos >= end) break;
                if (bed && bed_overlap(bed, sam_hdr_tid2name(h, last_tid), last_pos, last_pos + 1) == 0)
                    continue;
                fputs(sam_hdr_tid2name(h, last_tid), file_out);
                fprintf(file_out, "\t%"PRIhts_pos, last_pos+1);
                for (i = 0; i < n; i++)
                    fputc('\t', file_out), fputc('0', file_out);
                fputc('\n', file_out);
            }
            last_tid++;
            last_pos = -1;
            if (all < 2 || reg)
                break;
        }
    }

depth_end:
    if (fclose(file_out) != 0) {
        if (status == EXIT_SUCCESS) {
            print_error_errno("depth", "error on closing \"%s\"",
                              (output_file && strcmp(output_file, "-") != 0
                               ? output_file : "stdout"));
            status = EXIT_FAILURE;
        }
    }

    for (i = 0; i < n && data[i]; ++i) {
        sam_hdr_destroy(data[i]->hdr);
        if (data[i]->fp) sam_close(data[i]->fp);
        hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); free(reg);
    if (bed) bed_destroy(bed);
    if ( file_list )
    {
        for (i=0; i<n; i++) free(fn[i]);
        free(fn);
    }
    sam_global_args_free(&ga);
    return status;
}

#ifdef _MAIN_BAM2DEPTH
int main(int argc, char *argv[])
{
    return main_depth(argc, argv);
}
#endif
