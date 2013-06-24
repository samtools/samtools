#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "sam.h"

/*****************
 * callback APIs *
 *****************/

int bam_pileup_file(bamFile fp, int mask, bam_pileup_f func, void *func_data)
{
	bam_plbuf_t *buf;
	int ret;
	bam1_t *b;
	b = bam_init1();
	buf = bam_plbuf_init(func, func_data);
	bam_plbuf_set_mask(buf, mask);
	while ((ret = bam_read1(fp, b)) >= 0)
		bam_plbuf_push(b, buf);
	bam_plbuf_push(0, buf);
	bam_plbuf_destroy(buf);
	bam_destroy1(b);
	return 0;
}

void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask)
{
	bam_plp_set_mask(buf->iter, mask);
}

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
