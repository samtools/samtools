/*
    Author: Pierre Lindenbaum PhD @yokofakun
            Institut du Thorax - Nantes - France

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
#include <stdlib.h>
#include <stdio.h>
#include "sam_dynreadfilter.h"

/* The following dynamic read filter select reads containing an EcoRI site */

static const char* EcoRI ="GAATTC";
static const int LEN_ECORI = 6;
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN"; //copied from htslib/hts.c

int filterexample2_initialize(SamDynReadFilterPtr hook ,const bam_hdr_t *header) {
	return 1;
	}

int filterexample2_accept(SamDynReadFilterPtr hook ,const bam_hdr_t *header, bam1_t *b) {
        int i,j; 
        uint8_t *seq = bam_get_seq(b);
        if (b->core.l_qseq < 6) return 0;
        if (seq == NULL) return 0;
	for (i = 0; i +LEN_ECORI <= b->core.l_qseq; i++) {
               for(j=0;j< LEN_ECORI;j++) {
		    char c = seq_nt16_str[bam_seqi(seq,i+j)];
		    if (EcoRI[j] != c) break;
                    }
                if (j==LEN_ECORI) return 1;
		}
        return 0;
	}

int filterexample2_dispose(SamDynReadFilterPtr hook ,const bam_hdr_t *header, bam1_t *b) {
        return 1;
	}
