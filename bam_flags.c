/*  bam_flags.c -- flags subcommand.

    Copyright (C) 2013-2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
#include <htslib/sam.h>

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About: Convert between textual and numeric flag representation\n");
    fprintf(stderr, "Usage: samtools flags INT|STR[,...]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Flags:\n");
    fprintf(stderr, "\t0x%x\tPAIRED        .. paired-end (or multiple-segment) sequencing technology\n", BAM_FPAIRED);
    fprintf(stderr, "\t0x%x\tPROPER_PAIR   .. each segment properly aligned according to the aligner\n", BAM_FPROPER_PAIR);
    fprintf(stderr, "\t0x%x\tUNMAP         .. segment unmapped\n", BAM_FUNMAP);
    fprintf(stderr, "\t0x%x\tMUNMAP        .. next segment in the template unmapped\n", BAM_FMUNMAP);
    fprintf(stderr, "\t0x%x\tREVERSE       .. SEQ is reverse complemented\n", BAM_FREVERSE);
    fprintf(stderr, "\t0x%x\tMREVERSE      .. SEQ of the next segment in the template is reversed\n", BAM_FMREVERSE);
    fprintf(stderr, "\t0x%x\tREAD1         .. the first segment in the template\n", BAM_FREAD1);
    fprintf(stderr, "\t0x%x\tREAD2         .. the last segment in the template\n", BAM_FREAD2);
    fprintf(stderr, "\t0x%x\tSECONDARY     .. secondary alignment\n", BAM_FSECONDARY);
    fprintf(stderr, "\t0x%x\tQCFAIL        .. not passing quality controls\n", BAM_FQCFAIL);
    fprintf(stderr, "\t0x%x\tDUP           .. PCR or optical duplicate\n", BAM_FDUP);
    fprintf(stderr, "\t0x%x\tSUPPLEMENTARY .. supplementary alignment\n", BAM_FSUPPLEMENTARY);
    fprintf(stderr, "\n");
}


int main_flags(int argc, char *argv[])
{
    if ( argc!=2 ) usage();
    else
    {
        int mask = bam_str2flag(argv[1]);
        if ( mask<0 ) { fprintf(stderr,"Error: Could not parse \"%s\"\n", argv[1]); usage(); return 1; }
        printf("0x%x\t%d\t%s\n", mask, mask, bam_flag2str(mask));
    }
    return 0;
}

