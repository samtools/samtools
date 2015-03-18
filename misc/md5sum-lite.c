/*  md5sum-lite.c -- Basic md5sum implementation.

    Copyright (C) 2008, 2009 Genome Research Ltd.

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
#include <string.h>
#include "htslib/hts.h"

static void md5_one(const char *fn)
{
    unsigned char buf[4096], digest[16];
    char hex[33];
    hts_md5_context *md5;
    int l;
    FILE *fp;

    fp = strcmp(fn, "-")? fopen(fn, "r") : stdin;
    if (fp == 0) {
        fprintf(stderr, "md5sum: %s: No such file or directory\n", fn);
        exit(1);
    }
    if (!(md5 = hts_md5_init())) {
        fprintf(stderr, "md5sum: %s: Failed to allocate md5 buffer\n", fn);
        exit(1);
    }
    while ((l = fread(buf, 1, 4096, fp)) > 0)
        hts_md5_update(md5, buf, l);
    hts_md5_final(digest, md5);
    if (fp != stdin) fclose(fp);
    hts_md5_hex(hex, digest);
    printf("%s  %s\n", hex, fn);
    hts_md5_destroy(md5);
}
int main(int argc, char *argv[])
{
    int i;
    if (argc == 1) md5_one("-");
    else for (i = 1; i < argc; ++i) md5_one(argv[i]);
    return 0;
}
