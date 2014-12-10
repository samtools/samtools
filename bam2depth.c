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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "samtools.h"

typedef struct {     // auxiliary data structure
    samFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps

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

int main_depth(int argc, char *argv[])
{
    int i, n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, status = EXIT_SUCCESS, nfiles;
    const bam_pileup1_t **plp;
    char *reg = 0; // specified region
    void *bed = 0; // BED data structure
    char *file_list = NULL, **fn = NULL;
    bam_hdr_t *h = NULL; // BAM header of the 1st input
    aux_t **data;
    bam_mplp_t mplp;

    // parse the command line
    while ((n = getopt(argc, argv, "r:b:q:Q:l:f:")) >= 0) {
        switch (n) {
            case 'l': min_len = atoi(optarg); break; // minimum query length
            case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
            case 'b':
                bed = bed_read(optarg); // BED or position list file can be parsed now
                if (!bed) { print_error_errno("Could not read file \"%s\"", optarg); return 1; }
                break;
            case 'q': baseQ = atoi(optarg); break;   // base quality threshold
            case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
            case 'f': file_list = optarg; break;
        }
    }
    if (optind == argc && !file_list) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: samtools depth [options] in1.bam [in2.bam [...]]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -b <bed>            list of positions or regions\n");
        fprintf(stderr, "   -f <list>           list of input BAM filenames, one per line [null]\n");
        fprintf(stderr, "   -l <int>            read length threshold (ignore reads shorter than <int>)\n");
        fprintf(stderr, "   -q <int>            base quality threshold\n");
        fprintf(stderr, "   -Q <int>            mapping quality threshold\n");
        fprintf(stderr, "   -r <chr:from-to>    region\n");
        fprintf(stderr, "\n");
        return 1;
    }

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
    beg = 0; end = 1<<30;  // set the default region
    for (i = 0; i < n; ++i) {
        data[i] = calloc(1, sizeof(aux_t));
        data[i]->fp = sam_open(argv[optind+i], "r"); // open BAM
        if (data[i]->fp == NULL) {
            print_error_errno("Could not open \"%s\"", argv[optind+i]);
            status = EXIT_FAILURE;
            goto depth_end;
        }
        data[i]->min_mapQ = mapQ;                    // set the mapQ filter
        data[i]->min_len  = min_len;                 // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        if (reg) { // if a region is specified
            hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
            if (idx == NULL) {
                print_error("can't load index for \"%s\"", argv[optind+i]);
                status = EXIT_FAILURE;
                goto depth_end;
            }
            data[i]->iter = sam_itr_querys(idx, data[i]->hdr, reg); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
            if (data[i]->iter == NULL) {
                print_error("can't parse region \"%s\"", reg);
                status = EXIT_FAILURE;
                goto depth_end;
            }
        }
    }

    h = data[0]->hdr; // easy access to the header of the 1st BAM
    if (reg) {
        beg = data[0]->iter->beg; // and to the parsed region coordinates
        end = data[0]->iter->end;
    }

    // the core multi-pileup loop
    mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
    n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(n, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
        if (pos < beg || pos >= end) continue; // out of range; skip
        if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
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
    free(n_plp); free(plp);
    bam_mplp_destroy(mplp);

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
    return status;
}

#ifdef _MAIN_BAM2DEPTH
int main(int argc, char *argv[])
{
    return main_depth(argc, argv);
}
#endif
