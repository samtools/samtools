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
DEALINGS IN THE SOFTWARE.

History:

  * 2016-01-12: Pierre Lindenbaum @yokofakun : added options -o -n 

*/

#include <config.h>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <getopt.h>
#include <htslib/faidx.h>

#define DEFAULT_FASTA_LINE_LEN 60

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
        fprintf(stderr, "Usage:\n\n");
        fprintf(stderr, " Build Index:\n\n  samtools faidx <file.fa|file.fa.gz>\n\n");
        fprintf(stderr, " Query Index:\n\n  samtools faidx [-n|--length line_len number] [-o|--output output.fa] <file.fa|file.fa.gz> [<reg> [...]]\n\n",DEFAULT_FASTA_LINE_LEN);
        fprintf(stderr, "\n");
    }
    exit(-1);
}


int faidx_main(int argc, char *argv[])
{
    int c;
    int line_len = DEFAULT_FASTA_LINE_LEN ;/* fasta line len */
    char* output_file = NULL; /* output file (default is stdout ) */
    FILE* file_out = stdout;/* output stream */
    
    static const struct option lopts[] = {
        { "output", required_argument, NULL, 'o' },
        { "help",   no_argument,       NULL, 'h' },
        { "length", required_argument, NULL, 'n' },
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "ho:n:", lopts, NULL)) >= 0) {
        switch (c) {
            case 'o': output_file = optarg; break;
            case 'n': line_len = atoi(optarg);
                      if(line_len<1) {
                        fprintf(stderr,"[faidx] bad line length '%s', using default:%d\n",optarg,DEFAULT_FASTA_LINE_LEN);
                        line_len= DEFAULT_FASTA_LINE_LEN ;
                        }
                      break;
            case '?': fprintf(stderr,"[faidx] unknown option '%c'\n",optopt); return EXIT_FAILURE;break;
            case 'h': error(NULL);
            default:  break;
        }
    }

    
    if ( argc==optind )
        error(NULL);
    if ( optind+1 == argc )
    {
        if (fai_build(argv[optind]) != 0) {
            error("Could not build fai index %s.fai\n", argv[optind]);
        }
        return 0;
    }

    faidx_t *fai = fai_load(argv[optind]);
    if ( !fai ) error("[faidx] Could not load fai index of %s\n", argv[optind]);
    
    /** output file provided by user */
    if( output_file != NULL ) {
          if( strcmp( output_file, argv[optind] ) == 0 ) {
               fprintf(stderr,"[faidx] same input/output : %s\n", output_file);
               return EXIT_FAILURE;
               }
          file_out = fopen( output_file, "w" );
          if( file_out == NULL) {
               fprintf(stderr,"[faidx] Cannot open \"%s\" for writing :%s.\n", output_file, strerror(errno) );
               return EXIT_FAILURE;
               }
          }
    while ( ++optind<argc )
    {
        fprintf(file_out, ">%s\n", argv[optind]);
        int i, j, seq_len;
        char *seq = fai_fetch(fai, argv[optind], &seq_len);
        if ( seq_len < 0 ) error("[faidx] Failed to fetch sequence in %s\n", argv[optind]);
        for (i=0; i<seq_len; i+= line_len)
        {
            for (j=0; j<line_len && i+j<seq_len; j++)
                fputc(seq[i+j],file_out);
            fputc('\n',file_out);
        }
        free(seq);
    }
    fai_destroy(fai);
    fflush(file_out);
    if( output_file != NULL) fclose(file_out);
    return 0;
}

