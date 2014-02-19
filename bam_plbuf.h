/* based on code from bam.h */

#ifndef BAM_PLBUF_H
#define BAM_PLBUF_H

#include <htslib/sam.h>

#ifndef BAM_PILEUP_F_DEFINED
#define BAM_PILEUP_F_DEFINED
typedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
#endif //BAM_PILEUP_F_DEFINED

typedef struct {
	bam_plp_t iter;
	bam_pileup_f func;
	void *data;
} bam_plbuf_t;

#ifdef __cplusplus
extern "C" {
#endif
	void bam_plbuf_reset(bam_plbuf_t *buf);
	
	bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data);
	
	void bam_plbuf_destroy(bam_plbuf_t *buf);
	
	int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf);
#ifdef __cplusplus
}
#endif

#endif // BAM_PLBUF_H
