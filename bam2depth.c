/*  bam2depth.c -- depth subcommand.

    Copyright (C) 2011, 2012 Broad Institute.
    Copyright (C) 2012-2014 Genome Research Ltd.

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
#include "sam_opts.h"

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps
int bed_query(const void *_h, const char *chr, int pos, int *beg, int *end);

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1)
    {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
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
    fprintf(stderr, "   -f <list>           list of input BAM filenames, one per line [null]\n");
    fprintf(stderr, "   -l <int>            read length threshold (ignore reads shorter than <int>) [0]\n");
    fprintf(stderr, "   -d/-m <int>         maximum coverage depth [8000]\n");  // the htslib's default
    fprintf(stderr, "   -q <int>            base quality threshold [0]\n");
    fprintf(stderr, "   -Q <int>            mapping quality threshold [0]\n");
    fprintf(stderr, "   -r <chr:from-to>    region\n");

    sam_global_opt_help(stderr, "-.--.-");

    fprintf(stderr, "\n");
    fprintf(stderr, "The output is a simple tab-separated table with three columns: reference name,\n");
    fprintf(stderr, "position, and coverage depth.  Note that positions with zero coverage may be\n");
    fprintf(stderr, "omitted by default; see the -a option.\n");
    fprintf(stderr, "\n");

    return 1;
}

int main_depth(int argc, char *argv[])
{
    int i, n, tid, reg_tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = 0;
    int all = 0, status = EXIT_SUCCESS, nfiles, max_depth = -1;
    const bam_pileup1_t **plp;
    char *reg = 0; // specified region
    void *bed = 0; // BED data structure
    char *file_list = NULL, **fn = NULL;
    bam_hdr_t *h = NULL; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;
    int last_pos = -1, last_tid = -1, ret;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-'),
        { NULL, 0, NULL, 0 }
    };

    // parse the command line
    while ((n = getopt_long(argc, argv, "r:b:q:Q:l:f:am:d:", lopts, NULL)) >= 0) {
        switch (n) {
            case 'l': min_len = atoi(optarg); break; // minimum query length
            case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
            case 'b':
                bed = bed_read(optarg); // BED or position list file can be parsed now
                if (!bed) { print_error_errno("depth", "Could not read file \"%s\"", optarg); return 1; }
                break;
            case 'q': baseQ = atoi(optarg); break;   // base quality threshold
            case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
            case 'f': file_list = optarg; break;
            case 'a': all++; break;
            case 'd': case 'm': max_depth = atoi(optarg); break; // maximum coverage depth
            default:  if (parse_sam_global_opt(n, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
            case '?': return usage();
        }
    }
    if (optind == argc && !file_list)
        return usage();

    // initialize the auxiliary data structures
    if (file_list)
    {
        if ( read_file_list(file_list,&nfiles,&fn) ) return 1;
        n = nfiles;
        argv = fn;
        optind = 0;
    }
    else
        n = argc - optind; // the number of BAMs on the command line
    data = calloc(n, sizeof(aux_t*)); // data[i] for the i-th input
    reg_tid = 0; beg = 0; end = INT_MAX;  // set the default region
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
            fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
            return 1;
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            return 1;
        }
        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = min_len;                 // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        if (data[i]->hdr == NULL) {
            fprintf(stderr, "Couldn't read header for \"%s\"\n",
                    argv[optind+i]);
            status = EXIT_FAILURE;
            goto depth_end;
        }
        if (reg) { // if a region is specified
            hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
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
    n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        if (tid >= h->n_targets) continue;     // diff number of @SQ lines per file?
        if (all) {
            while (tid > last_tid) {
                if (last_tid >= 0 && !reg) {
                    // Deal with remainder or entirety of last tid.
                    while (++last_pos < h->target_len[last_tid]) {
                        // Horribly inefficient, but the bed API is an obfuscated black box.
                        if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                            continue;
                        fputs(h->target_name[last_tid], stdout); printf("\t%d", last_pos+1);
                        for (i = 0; i < n; i++)
                            putchar('\t'), putchar('0');
                        putchar('\n');
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
                if (bed && bed_overlap(bed, h->target_name[tid], last_pos, last_pos + 1) == 0)
                    continue;
                fputs(h->target_name[tid], stdout); printf("\t%d", last_pos+1);
                for (i = 0; i < n; i++)
                    putchar('\t'), putchar('0');
                putchar('\n');
            }

            last_tid = tid;
            last_pos = pos;
        }
        if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue;
        fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
        for (i = 0; i < n; ++i) { // base level filters have to go here
            int j, m = 0;
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
                if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
                else if (bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
            }
            printf("\t%d", n_plp[i] - m); // this the depth to output
        }
        putchar('\n');
    }
    if (ret < 0) status = EXIT_FAILURE;
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);

    if (all) {
        // Handle terminating region
        if (last_tid < 0 && reg && all > 1) {
            last_tid = reg_tid;
            last_pos = beg-1;
        }
        while (last_tid >= 0 && last_tid < h->n_targets) {
            while (++last_pos < h->target_len[last_tid]) {
                if (last_pos >= end) break;
                if (bed && bed_overlap(bed, h->target_name[last_tid], last_pos, last_pos + 1) == 0)
                    continue;
                fputs(h->target_name[last_tid], stdout); printf("\t%d", last_pos+1);
                for (i = 0; i < n; i++)
                    putchar('\t'), putchar('0');
                putchar('\n');
            }
            last_tid++;
            last_pos = -1;
            if (all < 2 || reg)
                break;
        }
    }

depth_end:
    for (i = 0; i < n && data[i]; ++i) {
        bam_hdr_destroy(data[i]->hdr);
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
