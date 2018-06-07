/*  fqidx.c -- fqidx subcommand.

    Copyright (C) 2008, 2009, 2013, 2016 Genome Research Ltd.
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

#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <htslib/fqidx.h>
#include "samtools.h"

static int usage(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools fqidx <file.fq|file.fq.gz> [<reg> [...]]\n");
    return exit_status;
}

int fqidx_main(int argc, char *argv[])
{
    int c;
    while((c  = getopt(argc, argv, "h")) >= 0)
    {
        switch(c)
        {
            case 'h':
                return usage(stdout, EXIT_SUCCESS);

            default:
                return usage(stderr, EXIT_FAILURE);
        }
    }
    if ( argc==optind )
        return usage(stdout, EXIT_SUCCESS);
    if ( argc==2 )
    {
        if (fqi_build(argv[optind]) != 0) {
            fprintf(stderr, "Could not build fqi index %s.fqi\n", argv[optind]);
            return EXIT_FAILURE;
        }
        return 0;
    }

    fqidx_t *fqi = fqi_load(argv[optind]);
    if ( !fqi ) {
        fprintf(stderr, "Could not load fqi index of %s\n", argv[optind]);
        return EXIT_FAILURE;
    }

    int exit_status = EXIT_SUCCESS;

    while ( ++optind<argc && exit_status == EXIT_SUCCESS)
    {
        printf("@%s\n", argv[optind]);
        int seq_len, qual_len;
        char *seq = fqi_fetch(fqi, argv[optind], &seq_len);
        if ( seq_len < 0 ) {
            fprintf(stderr, "Failed to fetch sequence in %s\n", argv[optind]);
            exit_status = EXIT_FAILURE;
            break;
        }
        size_t i, seq_sz = seq_len;
        for (i=0; i<seq_sz; i+=60)
        {
            size_t len = i + 60 < seq_sz ? 60 : seq_sz - i;
            if (fwrite(seq + i, 1, len, stdout) < len ||
                putchar('\n') == EOF) {
                print_error_errno("fqidx", "failed to write output");
                exit_status = EXIT_FAILURE;
                break;
            }
        }
        free(seq);
        printf("+\n");
        char *qual = fqi_fetchqual(fqi, argv[optind], &qual_len);
        if ( qual_len < 0 ) {
            fprintf(stderr, "Failed to fetch sequence in %s\n", argv[optind]);
            exit_status = EXIT_FAILURE;
            break;
        }
        seq_sz = qual_len;
        for (i=0; i<seq_sz; i+=60)
        {
            size_t len = i + 60 < seq_sz ? 60 : seq_sz - i;
            if (fwrite(qual + i, 1, len, stdout) < len ||
                putchar('\n') == EOF) {
                print_error_errno("fqidx", "failed to write output");
                exit_status = EXIT_FAILURE;
                break;
            }
        }

        free(qual);
    }
    fqi_destroy(fqi);

    if (fflush(stdout) == EOF) {
        print_error_errno("fqidx", "failed to flush output");
        exit_status = EXIT_FAILURE;
    }

    return exit_status;
}
