/*  test/vcf-miniview.c -- minimal BCF/VCF viewer, for use by test harness.

    Copyright (C) 2014 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/vcf.h>

void usage()
{
    fprintf(stderr,
"Usage: vcf-miniview [view] [-f] [FILE]\n"
"Options:\n"
"  -f  Filters out ## headers and various fields, to simplify file comparison\n");
    exit(EXIT_FAILURE);
}

void fail(const char *message)
{
    fprintf(stderr, "vcf-miniview: %s\n", message);
    exit(EXIT_FAILURE);
}

void erase(kstring_t* str, const char *tag)
{
    char *begin, *end;

    if ((begin = strstr(str->s+1, tag)) == NULL) return;

    end = begin;
    while (*end && *end != '\t' && *end != ';') end++;
    if (begin[-1] == ';') begin--;

    memmove(begin, end, str->l - (end - str->s) + 1);
    str->l -= end - begin;
}

int main(int argc, char **argv)
{
    int optind, filter = 0;
    htsFile *in;
    bcf_hdr_t *hdr;
    bcf1_t *rec;
    kstring_t str = { 0, 0, NULL };

    optind = 1;
    if (optind < argc && strcmp(argv[optind], "view") == 0) optind++;
    if (optind < argc && strcmp(argv[optind], "-f") == 0) filter = 1, optind++;
    if (argc == 1 || argc - optind > 1) usage();

    if ((in = hts_open((optind < argc)? argv[optind] : "-", "r")) == NULL)
        fail("can't open input file");

    if ((hdr = bcf_hdr_read(in)) == NULL)
        fail("can't read header");

    bcf_hdr_format(hdr, 0, &str);
    if (filter) {
        char *fixed = strstr(str.s, "\n#CHROM");
        printf("%s", fixed? fixed+1 : str.s);
    }
    else printf("%s", str.s);

    rec = bcf_init();
    while (bcf_read(in, hdr, rec) >= 0) {
        str.l = 0;
        vcf_format(hdr, rec, &str);
        if (filter) {
            erase(&str, "IMF=");
            erase(&str, "DP=");
            erase(&str, "IDV=");
            erase(&str, "IMP=");
            erase(&str, "IS=");
            erase(&str, "VDB=");
            erase(&str, "SGB=");
            erase(&str, "MQB=");
            erase(&str, "BQB=");
            erase(&str, "RPB=");
            erase(&str, "MQ0F=");
            erase(&str, "MQSB=");
        }
        printf("%s", str.s);
    }

    free(str.s);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(in);
    return EXIT_SUCCESS;
}
