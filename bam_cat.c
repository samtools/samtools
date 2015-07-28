/*  bam_cat.c -- efficiently concatenates bam files.

    Copyright (C) 2008-2009, 2011-2013 Genome Research Ltd.
    Modified SAMtools work copyright (C) 2010 Illumina, Inc.

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

/*
bam_cat can be used to concatenate BAM files. Under special
circumstances, it can be used as an alternative to 'samtools merge' to
concatenate multiple sorted files into a single sorted file. For this
to work each file must be sorted, and the sorted files must be given
as command line arguments in order such that the final read in file i
is less than or equal to the first read in file i+1.

This code is derived from the bam_reheader function in samtools 0.1.8
and modified to perform concatenation by Chris Saunders on behalf of
Illumina.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "htslib/bgzf.h"
#include "bam.h"

#define BUF_SIZE 0x10000

#define GZIPID1 31
#define GZIPID2 139

#define BGZF_EMPTY_BLOCK_SIZE 28


int bam_cat(int nfn, char * const *fn, const bam_header_t *h, const char* outbam)
{
    BGZF *fp, *in = NULL;
    uint8_t *buf = NULL;
    uint8_t ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int es=BGZF_EMPTY_BLOCK_SIZE;
    int i;

    fp = strcmp(outbam, "-")? bgzf_open(outbam, "w") : bgzf_fdopen(fileno(stdout), "w");
    if (fp == 0) {
        fprintf(stderr, "[%s] ERROR: fail to open output file '%s'.\n", __func__, outbam);
        return 1;
    }
    if (h) bam_header_write(fp, h);

    buf = (uint8_t*) malloc(BUF_SIZE);
    if (!buf) {
        fprintf(stderr, "[%s] Couldn't allocate buffer\n", __func__);
        goto fail;
    }
    for(i = 0; i < nfn; ++i){
        bam_header_t *old;
        int len,j;

        in = strcmp(fn[i], "-")? bgzf_open(fn[i], "r") : bgzf_fdopen(fileno(stdin), "r");
        if (in == 0) {
            fprintf(stderr, "[%s] ERROR: fail to open file '%s'.\n", __func__, fn[i]);
            goto fail;
        }
        if (in->is_write) return -1;

        old = bam_header_read(in);
        if (old == NULL) {
            fprintf(stderr, "[%s] ERROR: couldn't read header for '%s'.\n",
                    __func__, fn[i]);
            goto fail;
        }
        if (h == 0 && i == 0) bam_header_write(fp, old);

        if (in->block_offset < in->block_length) {
            if (bgzf_write(fp, in->uncompressed_block + in->block_offset, in->block_length - in->block_offset) < 0) goto write_fail;
            if (bgzf_flush(fp) != 0) goto write_fail;
        }

        j=0;
        while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0) {
            if(len<es){
                int diff=es-len;
                if(j==0) {
                    fprintf(stderr, "[%s] ERROR: truncated file?: '%s'.\n", __func__, fn[i]);
                    goto fail;
                }
                if (bgzf_raw_write(fp, ebuf, len) < 0) goto write_fail;
                
                memcpy(ebuf,ebuf+len,diff);
                memcpy(ebuf+diff,buf,len);
            } else {
                if(j!=0) {
                    if (bgzf_raw_write(fp, ebuf, es) < 0) goto write_fail;
                }
                len-= es;
                memcpy(ebuf,buf+len,es);
                if (bgzf_raw_write(fp, buf, len) < 0) goto write_fail;
            }
            j=1;
        }

        /* check final gzip block */
        {
            const uint8_t gzip1=ebuf[0];
            const uint8_t gzip2=ebuf[1];
            const uint32_t isize=*((uint32_t*)(ebuf+es-4));
            if(((gzip1!=GZIPID1) || (gzip2!=GZIPID2)) || (isize!=0)) {
                fprintf(stderr, "[%s] WARNING: Unexpected block structure in file '%s'.", __func__, fn[i]);
                fprintf(stderr, " Possible output corruption.\n");
                if (bgzf_raw_write(fp, ebuf, es) < 0) goto write_fail;
            }
        }
        bam_header_destroy(old);
        bgzf_close(in);
        in = NULL;
    }
    free(buf);
    if (bgzf_close(fp) < 0) {
        fprintf(stderr, "[%s] Error on closing '%s'.\n", __func__, outbam);
        return -1;
    }
    return 0;

 write_fail:
    fprintf(stderr, "[%s] Error writing to '%s'.\n", __func__, outbam);
 fail:
    if (in) bgzf_close(in);
    if (fp) bgzf_close(fp);
    free(buf);
    return -1;
}



int main_cat(int argc, char *argv[])
{
    bam_header_t *h = 0;
    char *outfn = 0;
    int c, ret;
    while ((c = getopt(argc, argv, "h:o:")) >= 0) {
        switch (c) {
            case 'h': {
                tamFile fph = sam_open(optarg);
                if (fph == 0) {
                    fprintf(stderr, "[%s] ERROR: fail to read the header from '%s'.\n", __func__, argv[1]);
                    return 1;
                }
                h = sam_header_read(fph);
                if (h == NULL) {
                    fprintf(stderr,
                            "[%s] ERROR: failed to read the header for '%s'.\n",
                            __func__, argv[1]);
                    return 1;
                }
                sam_close(fph);
                break;
            }
            case 'o': outfn = strdup(optarg); break;
        }
    }
    if (argc - optind < 2) {
        fprintf(stderr, "Usage: samtools cat [-h header.sam] [-o out.bam] <in1.bam> <in2.bam> [...]\n");
        return 1;
    }
    ret = bam_cat(argc - optind, argv + optind, h, outfn? outfn : "-");
    free(outfn);
    return ret;
}
