#include <stdio.h>
#include "sam.h"
static int fetch_func(const bam1_t *b, void *data)
{
	samfile_t *fp = (samfile_t*)data;
	uint32_t *cigar = bam1_cigar(b);
	const bam1_core_t *c = &b->core;
	int i, l;
	if (b->core.tid < 0) return 0;
	for (i = l = 0; i < c->n_cigar; ++i) {
		int op = cigar[i]&0xf;
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
			l += cigar[i]>>4;
	}
	printf("%s\t%d\t%d\t%s\t%d\t%c\n", fp->header->target_name[c->tid],
		   c->pos, c->pos + l, bam1_qname(b), c->qual, (c->flag&BAM_FREVERSE)? '-' : '+');
	return 0;
}
int main(int argc, char *argv[])
{
	samfile_t *fp;
	if (argc == 1) {
		fprintf(stderr, "Usage: bam2bed <in.bam> [region]\n");
		return 1;
	}
	if ((fp = samopen(argv[1], "rb", 0)) == 0) {
		fprintf(stderr, "bam2bed: Fail to open BAM file %s\n", argv[1]);
		return 1;
	}
	if (argc == 2) { /* if a region is not specified */
		bam1_t *b = bam_init1();
		while (samread(fp, b) >= 0) fetch_func(b, fp);
		bam_destroy1(b);
	} else {
		int ref, beg, end;
		bam_index_t *idx;
		if ((idx = bam_index_load(argv[1])) == 0) {
			fprintf(stderr, "bam2bed: BAM indexing file is not available.\n");
			return 1;
		}
		bam_parse_region(fp->header, argv[2], &ref, &beg, &end);
		if (ref < 0) {
			fprintf(stderr, "bam2bed: Invalid region %s\n", argv[2]);
			return 1;
		}
		bam_fetch(fp->x.bam, idx, ref, beg, end, fp, fetch_func);
		bam_index_destroy(idx);
	}
	samclose(fp);
	return 0;
}
