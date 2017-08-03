#ifndef BEDIDX_H
#define BEDIDX_H

#define LIDX_SHIFT 13
#define ALL 0
#define FILTERED 1

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
inline int bed_end(void *reg_hash);
void *bed_insert(void *reg_hash, char *reg, int beg, int end);
void *bed_filter(void *reg_hash, char *reg);
const char* bed_get(void *reg_hash, int index, int filter);

#endif
