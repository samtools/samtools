/*  faidx.c -- faidx subcommand.

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
#include <string.h>

#include <htslib/faidx.h>
#include "samtools.h"

int print_fasta(faidx_t *fai, char * region){
  int seq_len;
  int ex = EXIT_SUCCESS;
  char *seq = fai_fetch(fai, region, &seq_len);
  if ( seq_len < 0 ) {
      fprintf(stderr, "Failed to fetch sequence in %s\n", region);
      ex = EXIT_FAILURE;
      return ex;
  }
  printf(">%s\n", region);
  size_t i, seq_sz = seq_len;
  for (i=0; i<seq_sz; i+=60)
  {
      size_t len = i + 60 < seq_sz ? 60 : seq_sz - i;
      if (fwrite(seq + i, 1, len, stdout) < len ||
          putchar('\n') == EOF) {
          print_error_errno("faidx", "failed to write output");
          ex = EXIT_FAILURE;
          break;
      }
  }
  free(seq);
  return ex;
}

char * read_string_from_fp(FILE* fp, size_t size){
//The size is extended by the input with the value of the provisional
    char *str;
    int ch;
    size_t len = 0;
    str = realloc(NULL, sizeof(char)*size);//size is start size
    if(!str)return str;
    while(EOF!=(ch=fgetc(fp)) && ch != '\n'){
        str[len++]=ch;
        if(len==size){
            str = realloc(str, sizeof(char)*(size+=16000));
            if(!str)return str;
        }
    }
    str[len++]='\0';

    return realloc(str, sizeof(char)*len);
}


static int usage(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools faidx <file.fa|file.fa.gz> [<reg> [...]]\n");
    return exit_status;
}

int faidx_main(int argc, char *argv[])
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
        if (fai_build(argv[optind]) != 0) {
            fprintf(stderr, "Could not build fai index %s.fai\n", argv[optind]);
            return EXIT_FAILURE;
        }
        return 0;
    }

    faidx_t *fai = fai_load(argv[optind]);
    if ( !fai ) {
        fprintf(stderr, "Could not load fai index of %s\n", argv[optind]);
        return EXIT_FAILURE;
    }

    int exit_status = EXIT_SUCCESS;

    while ( ++optind<argc && exit_status == EXIT_SUCCESS)
    {
    /**
     * This checks if argv[optind] is a file containing regions. The format
     * needs to be the same, i.e. space delimited. This allows for a large
     * number of regions. A file or a region can be passed on the command line.
     * This maintains the UI and solves the "Argument too long" issue.
     */
      if( access( argv[optind], F_OK ) != -1 ) {
        FILE * fn = fopen(argv[optind], "r");
        char * regions  = read_string_from_fp(fn, 2000);
        fclose(fn);

        char* token = strtok(regions, " ");
        while (token && exit_status == EXIT_SUCCESS) {
          exit_status = print_fasta(fai, token);
          token = strtok(NULL, " ");
        }
      }
      else{
        exit_status = print_fasta(fai, argv[optind]);
      }
    }
    fai_destroy(fai);

    if (fflush(stdout) == EOF) {
        print_error_errno("faidx", "failed to flush output");
        exit_status = EXIT_FAILURE;
    }

    return exit_status;
}
