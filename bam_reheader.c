/*  bam_reheader.c -- reheader subcommand.

    Copyright (C) 2010 Broad Institute.
    Copyright (C) 2012, 2013 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include "htslib/bgzf.h"
#include "bam.h"

#define BUF_SIZE 0x10000

int bam_reheader(BGZF *in, const bam_header_t *h, int fd)
{
    BGZF *fp;
    bam_header_t *old;
    ssize_t len;
    uint8_t *buf;
    if (in->is_write) return -1;
    buf = malloc(BUF_SIZE);
    old = bam_header_read(in);
    fp = bgzf_fdopen(fd, "w");
    bam_header_write(fp, h);
    if (in->block_offset < in->block_length) {
        bgzf_write(fp, in->uncompressed_block + in->block_offset, in->block_length - in->block_offset);
        bgzf_flush(fp);
    }
    while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0)
        bgzf_raw_write(fp, buf, len);
    free(buf);
    fp->block_offset = in->block_offset = 0;
    bgzf_close(fp);
    return 0;
}

int main_reheader(int argc, char *argv[])
{
    bam_header_t *h;
    BGZF *in;
    if (argc != 3) {
        fprintf(stderr, "Usage: samtools reheader <in.header.sam> <in.bam>\n");
        return 1;
    }
    { // read the header
        tamFile fph = sam_open(argv[1]);
        if (fph == 0) {
            fprintf(stderr, "[%s] fail to read the header from %s.\n", __func__, argv[1]);
            return 1;
        }
        h = sam_header_read(fph);
        sam_close(fph);
    }
    in = strcmp(argv[2], "-")? bgzf_open(argv[2], "r") : bgzf_fdopen(fileno(stdin), "r");
    if (in == 0) {
        fprintf(stderr, "[%s] fail to open file %s.\n", __func__, argv[2]);
        return 1;
    }
    bam_reheader(in, h, fileno(stdout));
    bgzf_close(in);
    return 0;
}
