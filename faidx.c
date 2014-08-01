/*  faidx.c -- faidx subcommand.

    Copyright (C) 2008, 2009, 2013 Genome Research Ltd.
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

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <stdarg.h>
#include <htslib/faidx.h>

static void error(const char *format, ...)
{
    if ( format )
    {
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    else
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   samtools faidx <file.fa|file.fa.gz> [<reg> [...]]\n");
        fprintf(stderr, "\n");
    }
    exit(-1);
}


int faidx_main(int argc, char *argv[])
{
    int c;
    while((c  = getopt(argc, argv, "h")) >= 0)
    {
        switch(c)
        {
            case 'h':
            default:
                error(NULL);
        }
    }
    if ( argc==optind )
        error(NULL);
    if ( argc==2 )
    {
        fai_build(argv[optind]);
        return 0;
    }

    faidx_t *fai = fai_load(argv[optind]);
    if ( !fai ) error("Could not load fai index of %s\n", argv[optind]);

    while ( ++optind<argc )
    {
        printf(">%s\n", argv[optind]);
        int i, j, seq_len;
        char *seq = fai_fetch(fai, argv[optind], &seq_len);
        if ( seq_len < 0 ) error("Failed to fetch sequence in %s\n", argv[optind]);
        for (i=0; i<seq_len; i+=60)
        {
            for (j=0; j<60 && i+j<seq_len; j++)
                putchar(seq[i+j]);
            putchar('\n');
        }
        free(seq);
    }
    fai_destroy(fai);

    return 0;
}

