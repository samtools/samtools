#ifndef BAM_SAM_H
#define BAM_SAM_H

#include "bam.h"

typedef struct {
	int type;
	union {
		tamFile tamr;
		bamFile bam;
		FILE *tamw;
	} x;
	bam_header_t *header;
} samfile_t;

#ifdef __cplusplus
extern "C" {
#endif

	// mode can be: r/w/rb/wb. On writing, aux points to bam_header_t; on reading, aux points to the name of fn_list for SAM
	samfile_t *samopen(const char *fn, const char *mode, const void *aux);
	void samclose(samfile_t *fp);
	int samread(samfile_t *fp, bam1_t *b);
	int samwrite(samfile_t *fp, const bam1_t *b);
	int sampileup(samfile_t *fp, int mask, int min_mapQ, bam_pileup_f func, void *func_data);

#ifdef __cplusplus
}
#endif

#endif
