#include <unistd.h>
#include <assert.h>
#include "bam.h"

typedef struct {
	long long n_reads, n_mapped, n_pair_all, n_pair_map, n_pair_good;
	long long n_sgltn, n_read1, n_read2;
	long long n_qcfail, n_dup;
	long long n_diffchr, n_diffhigh;
} bam_flagstat_t;

#define flagstat_loop(s, c) do {										\
		++(s)->n_reads;													\
		if ((c)->flag & BAM_FPAIRED) {									\
			++(s)->n_pair_all;											\
			if ((c)->flag & BAM_FPROPER_PAIR) ++(s)->n_pair_good;		\
			if ((c)->flag & BAM_FREAD1) ++(s)->n_read1;					\
			if ((c)->flag & BAM_FREAD2) ++(s)->n_read2;					\
			if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP)) ++(s)->n_sgltn;	\
			if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP)) { \
				++(s)->n_pair_map;										\
				if ((c)->mtid != (c)->tid) {							\
					++(s)->n_diffchr;									\
					if ((c)->qual >= 5) ++(s)->n_diffhigh;				\
				}														\
			}															\
		}																\
		if (!((c)->flag & BAM_FUNMAP)) ++(s)->n_mapped;					\
		if ((c)->flag & BAM_FQCFAIL) ++(s)->n_qcfail;					\
		if ((c)->flag & BAM_FDUP) ++(s)->n_dup;							\
	} while (0)

bam_flagstat_t *bam_flagstat_core(bamFile fp)
{
	bam_flagstat_t *s;
	bam1_t *b;
	bam1_core_t *c;
	int ret;
	s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));
	b = bam_init1();
	c = &b->core;
	while ((ret = bam_read1(fp, b)) >= 0)
		flagstat_loop(s, c);
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");
	return s;
}
int bam_flagstat(int argc, char *argv[])
{
	bamFile fp;
	bam_header_t *header;
	bam_flagstat_t *s;
	if (argc == optind) {
		fprintf(stderr, "Usage: samtools flagstat <in.bam>\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);
	s = bam_flagstat_core(fp);
	printf("%lld in total\n", s->n_reads);
	printf("%lld QC failure\n", s->n_qcfail);
	printf("%lld duplicates\n", s->n_dup);
	printf("%lld mapped (%.2f%%)\n", s->n_mapped, (float)s->n_mapped / s->n_reads * 100.0);
	printf("%lld paired in sequencing\n", s->n_pair_all);
	printf("%lld read1\n", s->n_read1);
	printf("%lld read2\n", s->n_read2);
	printf("%lld properly paired (%.2f%%)\n", s->n_pair_good, (float)s->n_pair_good / s->n_pair_all * 100.0);
	printf("%lld with itself and mate mapped\n", s->n_pair_map);
	printf("%lld singletons (%.2f%%)\n", s->n_sgltn, (float)s->n_sgltn / s->n_pair_all * 100.0);
	printf("%lld with mate mapped to a different chr\n", s->n_diffchr);
	printf("%lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh);
	free(s);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}
