/*  md5fa.c -- MD5 checksum utility.

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

#include <config.h>

#include <stdio.h>
#include <zlib.h>
#include "htslib/kseq.h"
#include "htslib/hts.h"

KSEQ_INIT(gzFile, gzread)

static void md5_one(const char *fn)
{
    hts_md5_context *md5_one, *md5_all;
    int l, i, k;
    gzFile fp;
    kseq_t *seq;
    unsigned char unordered[16], digest[16];
    char hex[33];

    for (l = 0; l < 16; ++l) unordered[l] = 0;
    fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) {
        fprintf(stderr, "md5fa: %s: No such file or directory\n", fn);
        exit(1);
    }

    if (!(md5_all = hts_md5_init()))
        exit(1);

    if (!(md5_one = hts_md5_init())) {
        hts_md5_destroy(md5_all);
        exit(1);
    }

    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        for (i = k = 0; i < seq->seq.l; ++i) {
            if (islower(seq->seq.s[i])) seq->seq.s[k++] = toupper(seq->seq.s[i]);
            else if (isupper(seq->seq.s[i])) seq->seq.s[k++] = seq->seq.s[i];
        }
        hts_md5_reset(md5_one);
        hts_md5_update(md5_one, (unsigned char*)seq->seq.s, k);
        hts_md5_final(digest, md5_one);
        hts_md5_hex(hex, digest);
        for (l = 0; l < 16; ++l)
            unordered[l] ^= digest[l];
        printf("%s  %s  %s\n", hex, fn, seq->name.s);
        hts_md5_update(md5_all, (unsigned char*)seq->seq.s, k);
    }
    hts_md5_final(digest, md5_all);
    kseq_destroy(seq);

    hts_md5_hex(hex, digest);
    printf("%s  %s  >ordered\n", hex, fn);
    hts_md5_hex(hex, unordered);
    printf("%s  %s  >unordered\n", hex, fn);

    hts_md5_destroy(md5_all);
    hts_md5_destroy(md5_one);
}

int main(int argc, char *argv[])
{
    int i;
    if (argc == 1) md5_one("-");
    else for (i = 1; i < argc; ++i) md5_one(argv[i]);
    return 0;
}
