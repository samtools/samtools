/*  bam_flags.c -- flags subcommand.

    Copyright (C) 2013-2014, 2021 Genome Research Ltd.

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

#include <config.h>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <stdarg.h>
#include <htslib/sam.h>
#include "samtools.h"

static void usage(FILE *fp)
{
    static const struct { int bit; const char *desc; } *fl, flags[] = {
        { BAM_FPAIRED, "paired-end / multiple-segment sequencing technology" },
        { BAM_FPROPER_PAIR, "each segment properly aligned according to aligner" },
        { BAM_FUNMAP, "segment unmapped" },
        { BAM_FMUNMAP, "next segment in the template unmapped" },
        { BAM_FREVERSE, "SEQ is reverse complemented" },
        { BAM_FMREVERSE, "SEQ of next segment in template is rev.complemented" },
        { BAM_FREAD1, "the first segment in the template" },
        { BAM_FREAD2, "the last segment in the template" },
        { BAM_FSECONDARY, "secondary alignment" },
        { BAM_FQCFAIL, "not passing quality controls or other filters" },
        { BAM_FDUP, "PCR or optical duplicate" },
        { BAM_FSUPPLEMENTARY, "supplementary alignment" },
        { 0, NULL }
    };

    fprintf(fp,
"About: Convert between textual and numeric flag representation\n"
"Usage: samtools flags FLAGS...\n"
"\n"
"Each FLAGS argument is either an INT (in decimal/hexadecimal/octal) representing\n"
"a combination of the following numeric flag values, or a comma-separated string\n"
"NAME,...,NAME representing a combination of the following flag names:\n"
"\n");
    for (fl = flags; fl->desc; fl++) {
        char *name = bam_flag2str(fl->bit);
        fprintf(fp, "%#6x %5d  %-15s%s\n", fl->bit, fl->bit, name, fl->desc);
        free(name);
    }
}


int main_flags(int argc, char *argv[])
{
    if ( argc < 2 ) { usage(stdout); return 0; }

    int i;
    for (i = 1; i < argc; i++)
    {
        int mask = bam_str2flag(argv[i]);
        if ( mask<0 ) { print_error("flags", "Could not parse \"%s\"", argv[i]); usage(stderr); return 1; }
        char *str = bam_flag2str(mask);
        printf("0x%x\t%d\t%s\n", mask, mask, str);
        free(str);
    }
    return 0;
}

