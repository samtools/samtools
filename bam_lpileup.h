/* based on code from bam.h */

#ifndef BAM_LPILEUP_H
#define BAM_LPILEUP_H


#include <htslib/sam.h>

struct __bam_lplbuf_t;
typedef struct __bam_lplbuf_t bam_lplbuf_t;

#ifndef BAM_PILEUP_F_DEFINED
#define BAM_PILEUP_F_DEFINED
typedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
#endif //BAM_PILEUP_F_DEFINED


#ifdef __cplusplus
extern "C" {
#endif
	void bam_lplbuf_reset(bam_lplbuf_t *buf);

	/*! @abstract  bam_plbuf_init() equivalent with level calculated. */
	bam_lplbuf_t *bam_lplbuf_init(bam_pileup_f func, void *data);

	/*! @abstract  bam_plbuf_destroy() equivalent with level calculated. */
	void bam_lplbuf_destroy(bam_lplbuf_t *tv);

	/*! @abstract  bam_plbuf_push() equivalent with level calculated. */
	int bam_lplbuf_push(const bam1_t *b, bam_lplbuf_t *buf);
#ifdef __cplusplus
}
#endif

#endif // BAM_LPILEUP_H
