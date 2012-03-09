#include <stdio.h>
#include "sam.h"

typedef struct {
	int beg, end;
	samfile_t *in;
} tmpstruct_t;

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
	bam_plbuf_t *buf = (bam_plbuf_t*)data;
	bam_plbuf_push(b, buf);
	return 0;
}
// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	tmpstruct_t *tmp = (tmpstruct_t*)data;
	if ((int)pos >= tmp->beg && (int)pos < tmp->end)
		printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);
	return 0;
}

int main(int argc, char *argv[])
{
	tmpstruct_t tmp;
	if (argc == 1) {
		fprintf(stderr, "Usage: calDepth <in.bam> [region]\n");
		return 1;
	}
	tmp.beg = 0; tmp.end = 0x7fffffff;
	tmp.in = samopen(argv[1], "rb", 0);
	if (tmp.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);
		return 1;
	}
	if (argc == 2) { // if a region is not specified
		sampileup(tmp.in, -1, pileup_func, &tmp);
	} else {
		int ref;
		bam_index_t *idx;
		bam_plbuf_t *buf;
		idx = bam_index_load(argv[1]); // load BAM index
		if (idx == 0) {
			fprintf(stderr, "BAM indexing file is not available.\n");
			return 1;
		}
		bam_parse_region(tmp.in->header, argv[2], &ref, &tmp.beg, &tmp.end); // parse the region
		if (ref < 0) {
			fprintf(stderr, "Invalid region %s\n", argv[2]);
			return 1;
		}
		buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
		bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
	}
	samclose(tmp.in);
	return 0;
}
