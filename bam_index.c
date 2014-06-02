#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include <stdlib.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#define BAM_LIDX_SHIFT    14

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
	else bam_index_build(argv[1], BAM_LIDX_SHIFT);
	return 0;
}

int bam_idxstats(int argc, char *argv[])
{
	hts_idx_t* idx;
	bam_hdr_t* header;
	samFile* fp;

	if (argc < 2) {
		fprintf(stderr, "Usage: samtools idxstats <in.bam>\n");
		return 1;
	}
	fp = sam_open(argv[1], "r");
	if (fp == NULL) { fprintf(stderr, "[%s] fail to open BAM.\n", __func__); return 1; }
	header = sam_hdr_read(fp);
	sam_close(fp);
	idx = sam_index_load(fp, argv[1]);
	if (idx == NULL) { fprintf(stderr, "[%s] fail to load the index.\n", __func__); return 1; }

	int i;
	for (i = 0; i < header->n_targets; ++i) {
		// Print out contig name and length
		printf("%s\t%d", header->target_name[i], header->target_len[i]);
		// Now fetch info about it from the meta bin
		uint64_t u, v;
		hts_idx_get_stat(idx, i, &u, &v);
		printf("\t%" PRIu64 "\t%" PRIu64 "\n", u, v);
	}
	// Dump information about unmapped reads
	printf("*\t0\t0\t%" PRIu64 "\n", hts_idx_get_n_no_coor(idx));
	bam_hdr_destroy(header);
	hts_idx_destroy(idx);
	return 0;
}
