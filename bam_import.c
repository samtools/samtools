/*  bam_import.c -- SAM format parsing.

    Copyright (C) 2008-2013 Genome Research Ltd.

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

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "bam.h"
#include "htslib/kseq.h"

KSTREAM_INIT(gzFile, gzread, 16384)

bam_header_t *sam_header_read2(const char *fn)
{
    bam_header_t *header;
    int c, dret, n_targets = 0;
    gzFile fp;
    kstream_t *ks;
    kstring_t *str;
    kstring_t samstr = { 0, 0, NULL };
    if (fn == 0) return 0;
    fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
    if (fp == 0) return 0;
    ks = ks_init(fp);
    str = (kstring_t*)calloc(1, sizeof(kstring_t));
    while (ks_getuntil(ks, 0, str, &dret) > 0) {
        ksprintf(&samstr, "@SQ\tSN:%s", str->s);
        ks_getuntil(ks, 0, str, &dret);
        ksprintf(&samstr, "\tLN:%d\n", atoi(str->s));
        n_targets++;
        if (dret != '\n')
            while ((c = ks_getc(ks)) != '\n' && c != -1);
    }
    ks_destroy(ks);
    gzclose(fp);
    free(str->s); free(str);
    header = sam_hdr_parse(samstr.l, samstr.s? samstr.s : "");
    free(samstr.s);
    fprintf(stderr, "[sam_header_read2] %d sequences loaded.\n", n_targets);
    return header;
}
