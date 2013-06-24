#include "bam.h"

int bam_index_build2(const char *fn, const char *_fnidx)
{
	fprintf(stderr, "Samtools-htslib-API: bam_index_build2() not yet implemented\n");
	abort();
}

int bam_index(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: samtools index <in.bam> [out.index]\n");
		return 1;
	}
	if (argc >= 3) bam_index_build2(argv[1], argv[2]);
	else bam_index_build(argv[1]);
	return 0;
}

int bam_idxstats(int argc, char *argv[])
{
#if 0
	bam_index_t *idx;
	bam_header_t *header;
	bamFile fp;
	int i;
	if (argc < 2) {
		fprintf(stderr, "Usage: samtools idxstats <in.bam>\n");
		return 1;
	}
	fp = bam_open(argv[1], "r");
	if (fp == 0) { fprintf(stderr, "[%s] fail to open BAM.\n", __func__); return 1; }
	header = bam_header_read(fp);
	bam_close(fp);
	idx = bam_index_load(argv[1]);
	if (idx == 0) { fprintf(stderr, "[%s] fail to load the index.\n", __func__); return 1; }
	for (i = 0; i < idx->n; ++i) {
		khint_t k;
		khash_t(i) *h = idx->index[i];
		printf("%s\t%d", header->target_name[i], header->target_len[i]);
		k = kh_get(i, h, BAM_MAX_BIN);
		if (k != kh_end(h))
			printf("\t%llu\t%llu", (long long)kh_val(h, k).list[1].u, (long long)kh_val(h, k).list[1].v);
		else printf("\t0\t0");
		putchar('\n');
	}
	printf("*\t0\t0\t%llu\n", (long long)idx->n_no_coor);
	bam_header_destroy(header);
	bam_index_destroy(idx);
	return 0;
#else
	fprintf(stderr, "Samtools-htslib: bam_idxstats() not yet implemented\n");
	abort();
#endif
}

int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
{
	int ret;
	bam_iter_t iter;
	bam1_t *b;
	b = bam_init1();
	iter = bam_iter_query(idx, tid, beg, end);
	while ((ret = bam_iter_read(fp, iter, b)) >= 0) func(b, data);
	bam_iter_destroy(iter);
	bam_destroy1(b);
	return (ret == -1)? 0 : ret;
}
