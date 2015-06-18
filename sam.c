/*  sam.c -- format-neutral SAM/BAM API.

    Copyright (C) 2009, 2012-2015 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.

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

#include <string.h>
#include <unistd.h>
#include "htslib/faidx.h"
#include "sam.h"

int samthreads(samfile_t *fp, int n_threads, int n_sub_blks)
{
    if (hts_get_format(fp->file)->format != bam || !fp->is_write) return -1;
    bgzf_mt(fp->x.bam, n_threads, n_sub_blks);
    return 0;
}

samfile_t *samopen(const char *fn, const char *mode, const void *aux)
{
    // hts_open() is really sam_open(), except for #define games
    samFile *hts_fp = hts_open(fn, mode);
    if (hts_fp == NULL)  return NULL;

    samfile_t *fp = malloc(sizeof (samfile_t));
    fp->file = hts_fp;
    fp->x.bam = hts_fp->fp.bgzf;
    if (strchr(mode, 'r')) {
        if (aux) {
            if (hts_set_fai_filename(fp->file, aux) != 0) {
                sam_close(hts_fp);
                free(fp);
                return NULL;
            }
        }
        fp->header = sam_hdr_read(fp->file);  // samclose() will free this
        if (fp->header == NULL) {
            sam_close(hts_fp);
            free(fp);
            return NULL;
        }
        fp->is_write = 0;
        if (fp->header->n_targets == 0 && bam_verbose >= 1)
            fprintf(stderr, "[samopen] no @SQ lines in the header.\n");
    }
    else {
        enum htsExactFormat fmt = hts_get_format(fp->file)->format;
        fp->header = (bam_hdr_t *)aux;  // For writing, we won't free it
        fp->is_write = 1;
        if (!(fmt == text_format || fmt == sam) || strchr(mode, 'h')) sam_hdr_write(fp->file, fp->header);
    }

    return fp;
}

void samclose(samfile_t *fp)
{
    if (fp) {
        if (!fp->is_write && fp->header) bam_hdr_destroy(fp->header);
        sam_close(fp->file);
        free(fp);
    }
}

int samfetch(samfile_t *fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
{
    bam1_t *b = bam_init1();
    hts_itr_t *iter = sam_itr_queryi(idx, tid, beg, end);
    int ret;
    while ((ret = sam_itr_next(fp->file, iter, b)) >= 0) func(b, data);
    hts_itr_destroy(iter);
    bam_destroy1(b);
    return (ret == -1)? 0 : ret;
}

int sampileup(samfile_t *fp, int mask, bam_pileup_f func, void *func_data)
{
    bam_plbuf_t *buf;
    int ret;
    bam1_t *b;
    b = bam_init1();
    buf = bam_plbuf_init(func, func_data);
    if (mask < 0) mask = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    else mask |= BAM_FUNMAP;
    while ((ret = samread(fp, b)) >= 0) {
        // bam_plp_push() itself now filters out unmapped reads only
        if (b->core.flag & mask) b->core.flag |= BAM_FUNMAP;
        bam_plbuf_push(b, buf);
    }
    bam_plbuf_push(0, buf);
    bam_plbuf_destroy(buf);
    bam_destroy1(b);
    return 0;
}

char *samfaipath(const char *fn_ref)
{
    char *fn_list = 0;
    if (fn_ref == 0) return 0;
    fn_list = calloc(strlen(fn_ref) + 5, 1);
    strcat(strcpy(fn_list, fn_ref), ".fai");
    if (access(fn_list, R_OK) == -1) { // fn_list is unreadable
        if (access(fn_ref, R_OK) == -1) {
            fprintf(stderr, "[samfaipath] fail to read file %s.\n", fn_ref);
        } else {
            if (bam_verbose >= 3) fprintf(stderr, "[samfaipath] build FASTA index...\n");
            if (fai_build(fn_ref) == -1) {
                fprintf(stderr, "[samfaipath] fail to build FASTA index.\n");
                free(fn_list); fn_list = 0;
            }
        }
    }
    return fn_list;
}
