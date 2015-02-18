/*  test/test.c -- test harness utility routines.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Martin O. Pollard <mp15@sanger.ac.uk>

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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

#include "test.h"

void xfreopen(const char *path, const char *mode, FILE *stream)
{
    if (freopen(path, mode, stream) == NULL) {
        fprintf(stderr, __FILE__": error reopening %s: %s\n",
                path, strerror(errno));
        exit(2);
    }
}

void dump_hdr(const bam_hdr_t* hdr)
{
    printf("n_targets: %d\n", hdr->n_targets);
    printf("ignore_sam_err: %d\n", hdr->ignore_sam_err);
    printf("l_text: %u\n", hdr->l_text);
    printf("idx\ttarget_len\ttarget_name:\n");
    int32_t target;
    for (target = 0; target < hdr->n_targets; ++target) {
        printf("%d\t%u\t\"%s\"\n", target, hdr->target_len[target], hdr->target_name[target]);
    }
    printf("->text:(");
    fwrite((void*)hdr->text, (size_t) hdr->l_text, 1, stdout);
    printf(")\n");
}

void dump_read(const bam1_t* b) {
    printf("->core.tid:(%d)\n", b->core.tid);
    printf("->core.pos:(%d)\n", b->core.pos);
    printf("->core.bin:(%d)\n", b->core.bin);
    printf("->core.qual:(%d)\n", b->core.qual);
    printf("->core.l_qname:(%d)\n", b->core.l_qname);
    printf("->core.flag:(%d)\n", b->core.flag);
    printf("->core.n_cigar:(%d)\n", b->core.n_cigar);
    printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
    printf("->core.mtid:(%d)\n", b->core.mtid);
    printf("->core.mpos:(%d)\n", b->core.mpos);
    printf("->core.isize:(%d)\n", b->core.isize);
    if (b->data) {
        printf("->data:");
        int i;
        for (i = 0; i < b->l_data; ++i) {
            printf("%x ", b->data[i]);
        }
        printf("\n");
    }
    if (b->core.l_qname) {
        printf("qname: %s\n",bam_get_qname(b));
    }
    if (b->core.l_qseq) {
        printf("qseq:");
        int i;
        for (i = 0; i < b->core.l_qseq; ++i) {
            printf("%c",seq_nt16_str[seq_nt16_table[bam_seqi(bam_get_seq(b),i)]]);
        }
        printf("\n");
        printf("qual:");
        for (i = 0; i < b->core.l_qseq; ++i) {
            printf("%c",bam_get_qual(b)[i]);
        }
        printf("\n");
        
    }
    
    if (bam_get_l_aux(b)) {
        int i = 0;
        uint8_t* aux = bam_get_aux(b);
        
        while (i < bam_get_l_aux(b)) {
            printf("%.2s:%c:",aux+i,*(aux+i+2));
            i += 2;
            switch (*(aux+i)) {
                case 'Z':
                while (*(aux+1+i) != '\0') { putc(*(aux+1+i), stdout); ++i; }
                break;
            }
            putc('\n',stdout);
            ++i;++i;
        }
    }
    printf("\n");
}