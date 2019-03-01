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

/* The following dynamic read filter select clipped reads. $MIN_CIGAR_LEN can be used to set a minimal read length */


/* private data used by this filter */
typedef struct counter_t {
	long total; // number of reads seen
	long accepted; // number of accepted reads 
        int min_cigar_len;// min clip-length
	}Counter;


int filterexample1_initialize(SamDynReadFilterPtr hook ,const bam_hdr_t *header) {
	const char* env = getenv("MIN_CIGAR_LEN");
	Counter* counter =  (Counter*)calloc(1,sizeof(Counter));
        if (counter==NULL) return 0;
	hook->data = counter;
        if (env != NULL) counter->min_cigar_len = atoi(env);
	if (counter->min_cigar_len<1) counter->min_cigar_len = 1;
	return 1;
	}

int filterexample1_accept(SamDynReadFilterPtr hook ,const bam_hdr_t *header, bam1_t *b) {
        uint32_t *cigar;
        int32_t i;
	Counter* counter = (Counter*)hook->data;
	counter->total++;
        /* read must be mapped and number of cigar elements > 1*/
        if ((b)->core.flag&BAM_FUNMAP) return 0;
        if ((b)->core.n_cigar < 2 ) return 0;

        cigar = bam_get_cigar(b);
        /* no cigar ?*/
	if (cigar == NULL) return 0;
 
        /** loop over the cigar and search for a clipped cigar element with length>= min_len */
        for (i = 0; i< b->core.n_cigar ; ++i) {
            char c = bam_cigar_opchr(cigar[i]);
            if ( (c == 'S' || c == 'H' ) &&  bam_cigar_oplen(cigar[i]) >= counter->min_cigar_len ) {
                 break;
                 }
            }

        /* no clip found */
        if (i == b->core.n_cigar ) return 0;

	counter->accepted++;
        return 1;
	}

int filterexample1_dispose(SamDynReadFilterPtr hook ,const bam_hdr_t *header, bam1_t *b) {
        if (hook->data != NULL) {
                Counter* counter = (Counter*)hook->data;
                fprintf(stderr,"[FILTER EXAMPLE] read size = %d\n",counter->min_cigar_len);
                fprintf(stderr,"[FILTER EXAMPLE] read count = %ld\n",counter->total);
                fprintf(stderr,"[FILTER EXAMPLE] read accepted = %ld\n",counter->accepted);
		free(hook->data);
		hook->data = NULL;
		}
        return 1;
	}
