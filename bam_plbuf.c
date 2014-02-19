/* based on code from bam_pileup.c */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include "bam_plbuf.h"

/*****************
 * callback APIs *
 *****************/

void bam_plbuf_reset(bam_plbuf_t *buf)
{
	bam_plp_reset(buf->iter);
}

bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)
{
	bam_plbuf_t *buf;
	buf = calloc(1, sizeof(bam_plbuf_t));
	buf->iter = bam_plp_init(0, 0);
	buf->func = func;
	buf->data = data;
	return buf;
}

void bam_plbuf_destroy(bam_plbuf_t *buf)
{
	bam_plp_destroy(buf->iter);
	free(buf);
}

int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf)
{
	int ret, n_plp, tid, pos;
	const bam_pileup1_t *plp;
	ret = bam_plp_push(buf->iter, b);
	if (ret < 0) return ret;
	while ((plp = bam_plp_next(buf->iter, &tid, &pos, &n_plp)) != 0)
		buf->func(tid, pos, n_plp, plp, buf->data);
	return 0;
}
