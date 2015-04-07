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

#include <stdio.h>
#include <zlib.h>
#include "md5.h"
#include "htslib/kseq.h"

#define HEX_STR "0123456789abcdef"

KSEQ_INIT(gzFile, gzread)

static void md5_one(const char *fn)
{
    MD5_CTX md5_one, md5_all;
    int l, i, k;
    gzFile fp;
    kseq_t *seq;
    unsigned char unordered[16], digest[16];

    for (l = 0; l < 16; ++l) unordered[l] = 0;
    fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) {
        fprintf(stderr, "md5fa: %s: No such file or directory\n", fn);
        exit(1);
    }

    MD5Init(&md5_all);
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
        for (i = k = 0; i < seq->seq.l; ++i) {
            if (islower(seq->seq.s[i])) seq->seq.s[k++] = toupper(seq->seq.s[i]);
            else if (isupper(seq->seq.s[i])) seq->seq.s[k++] = seq->seq.s[i];
        }
        MD5Init(&md5_one);
        MD5Update(&md5_one, (unsigned char*)seq->seq.s, k);
        MD5Final(digest, &md5_one);
        for (l = 0; l < 16; ++l) {
            printf("%c%c", HEX_STR[digest[l]>>4&0xf], HEX_STR[digest[l]&0xf]);
            unordered[l] ^= digest[l];
        }
        printf("  %s  %s\n", fn, seq->name.s);
        MD5Update(&md5_all, (unsigned char*)seq->seq.s, k);
    }
    MD5Final(digest, &md5_all);
    kseq_destroy(seq);
    for (l = 0; l < 16; ++l)
        printf("%c%c", HEX_STR[digest[l]>>4&0xf], HEX_STR[digest[l]&0xf]);
    printf("  %s  >ordered\n", fn);
    for (l = 0; l < 16; ++l)
        printf("%c%c", HEX_STR[unordered[l]>>4&0xf], HEX_STR[unordered[l]&0xf]);
    printf("  %s  >unordered\n", fn);
}

int main(int argc, char *argv[])
{
    int i;
    if (argc == 1) md5_one("-");
    else for (i = 1; i < argc; ++i) md5_one(argv[i]);
    return 0;
}
