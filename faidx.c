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

