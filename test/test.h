#ifndef TEST_H
#define TEST_H

#include <stdio.h>
#include <htslib/sam.h>

void xfreopen(const char *path, const char *mode, FILE *stream);

void dump_hdr(const bam_hdr_t* hdr);
void dump_read(const bam1_t* read);

#endif
