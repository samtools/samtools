/*  faidx.c -- faidx subcommand.

    Copyright (C) 2008, 2009, 2013, 2016, 2018 Genome Research Ltd.
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
DEALINGS IN THE SOFTWARE.

History:

  * 2016-01-12: Pierre Lindenbaum @yokofakun : added options -o -n

*/

#include <config.h>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include "samtools.h"

#define DEFAULT_FASTA_LINE_LEN 60

static int usage(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools faidx <file.fa|file.fa.gz> [<reg> [...]]\n");
    fprintf(fp, "Option: \n"
                " -o, --output   FILE Write FASTA to file.\n"
                " -n, --length   INT  Length of FASTA sequence line. [60]\n"
                " -c, --continue      Continue after trying to retrieve missing region.\n"
                " -h, --help          This message.\n");
    return exit_status;
}

int faidx_main(int argc, char *argv[])
{
    int c, ignore_error = 0;
    int line_len = DEFAULT_FASTA_LINE_LEN ;/* fasta line len */
    char* output_file = NULL; /* output file (default is stdout ) */
    FILE* file_out = stdout;/* output stream */

    static const struct option lopts[] = {
        { "output", required_argument, NULL, 'o' },
        { "help",   no_argument,       NULL, 'h' },
        { "length", required_argument, NULL, 'n' },
        { "continue", no_argument,     NULL, 'c' },
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "ho:n:c", lopts, NULL)) >= 0) {
        switch (c) {
            case 'o': output_file = optarg; break;
            case 'n': line_len = atoi(optarg);
                      if(line_len<1) {
                        fprintf(stderr,"[faidx] bad line length '%s', using default:%d\n",optarg,DEFAULT_FASTA_LINE_LEN);
                        line_len= DEFAULT_FASTA_LINE_LEN ;
                        }
                      break;
            case 'c': ignore_error = 1; break;
            case '?': return usage(stderr, EXIT_FAILURE);
            case 'h': return usage(stdout, EXIT_SUCCESS);
            default:  break;
        }
    }

    if ( argc==optind )
        return usage(stdout, EXIT_SUCCESS);

    if ( optind+1 == argc )
    {
        if (fai_build(argv[optind]) != 0) {
            fprintf(stderr, "[faidx] Could not build fai index %s.fai\n", argv[optind]);
            return EXIT_FAILURE;
        }
        return 0;
    }

    faidx_t *fai = fai_load(argv[optind]);

    if ( !fai ) {
        fprintf(stderr, "[faidx] Could not load fai index of %s\n", argv[optind]);
        return EXIT_FAILURE;
    }

    /** output file provided by user */
    if( output_file != NULL ) {
        if( strcmp( output_file, argv[optind] ) == 0 ) {
            fprintf(stderr,"[faidx] Same input/output : %s\n", output_file);
            return EXIT_FAILURE;
        }

        file_out = fopen( output_file, "w" );

        if( file_out == NULL) {
            fprintf(stderr,"[faidx] Cannot open \"%s\" for writing :%s.\n", output_file, strerror(errno) );
            return EXIT_FAILURE;
        }
    }

    int exit_status = EXIT_SUCCESS;

    while ( ++optind<argc && exit_status == EXIT_SUCCESS)
    {
        printf(">%s\n", argv[optind]);
        int seq_len, beg, end;
        char *seq = fai_fetch(fai, argv[optind], &seq_len);

        if ( seq_len < 0 ) {
            fprintf(stderr, "[faidx] Failed to fetch sequence in %s\n", argv[optind]);

            if (ignore_error && seq_len == -2) {
                continue;
            } else {
                exit_status = EXIT_FAILURE;
                break;
            }
        } else if (seq_len == 0) {
            fprintf(stderr, "[faidx] Zero length sequence: %s\n", argv[optind]);
        } else if (hts_parse_reg(argv[optind], &beg, &end) && (end < INT_MAX) && (seq_len != end - beg)) {
            fprintf(stderr, "[faidx] Truncated sequence: %s\n", argv[optind]);
        }

        size_t i, seq_sz = seq_len;
        for (i=0; i<seq_sz; i+=line_len)
        {
            size_t len = i + line_len < seq_sz ? line_len : seq_sz - i;
            if (fwrite(seq + i, 1, len, file_out) < len ||
                fputc('\n', file_out) == EOF) {
                print_error_errno("faidx", "failed to write output");
                exit_status = EXIT_FAILURE;
                break;
            }
        }
        free(seq);
    }
    fai_destroy(fai);

    if (fflush(file_out) == EOF) {
        print_error_errno("faidx", "failed to flush output");
        exit_status = EXIT_FAILURE;
    }

    if( output_file != NULL) fclose(file_out);

    return exit_status;
}
