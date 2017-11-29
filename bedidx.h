#ifndef BEDIDX_H
#define BEDIDX_H

#include "htslib/hts.h"

#define LIDX_SHIFT 13
#define ALL 0
#define FILTERED 1

#define MIN(A,B) ( ( (A) < (B) ) ? (A) : (B) )
#define MAX(A,B) ( ( (A) > (B) ) ? (A) : (B) )

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void *bed_hash_regions(void *reg_hash, char **regs, int first, int last, int *op);
const char* bed_get(void *reg_hash, int index, int filter);
hts_reglist_t *bed_reglist(void *reg_hash, int filter, int *count_regs);

#endif
