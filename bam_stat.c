#include <unistd.h>
#include <assert.h>
#include "bam.h"

typedef struct {
	long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
	long long n_sgltn[2], n_read1[2], n_read2[2];
	long long n_dup[2];
	long long n_diffchr[2], n_diffhigh[2];
} bam_flagstat_t;

#define flagstat_loop(s, c) do {										\
		int w = ((c)->flag & BAM_FQCFAIL)? 1 : 0;	 					\
		++(s)->n_reads[w];												\
		if ((c)->flag & BAM_FPAIRED) {									\
			++(s)->n_pair_all[w];										\
			if ((c)->flag & BAM_FPROPER_PAIR) ++(s)->n_pair_good[w];	\
			if ((c)->flag & BAM_FREAD1) ++(s)->n_read1[w];				\
			if ((c)->flag & BAM_FREAD2) ++(s)->n_read2[w];				\
			if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP)) ++(s)->n_sgltn[w];	\
			if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP)) { \
				++(s)->n_pair_map[w];									\
				if ((c)->mtid != (c)->tid) {							\
					++(s)->n_diffchr[w];								\
					if ((c)->qual >= 5) ++(s)->n_diffhigh[w];			\
				}														\
			}															\
		}																\
		if (!((c)->flag & BAM_FUNMAP)) ++(s)->n_mapped[w];				\
		if ((c)->flag & BAM_FDUP) ++(s)->n_dup[w];						\
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
	printf("%lld + %lld in total (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
	printf("%lld + %lld duplicates\n", s->n_dup[0], s->n_dup[1]);
	printf("%lld + %lld mapped (%.2f%%:%.2f%%)\n", s->n_mapped[0], s->n_mapped[1], (float)s->n_mapped[0] / s->n_reads[0] * 100.0, (float)s->n_mapped[1] / s->n_reads[1] * 100.0);
	printf("%lld + %lld paired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
	printf("%lld + %lld read1\n", s->n_read1[0], s->n_read1[1]);
	printf("%lld + %lld read2\n", s->n_read2[0], s->n_read2[1]);
	printf("%lld + %lld properly paired (%.2f%%:%.2f%%)\n", s->n_pair_good[0], s->n_pair_good[1], (float)s->n_pair_good[0] / s->n_pair_all[0] * 100.0, (float)s->n_pair_good[1] / s->n_pair_all[1] * 100.0);
	printf("%lld + %lld with itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
	printf("%lld + %lld singletons (%.2f%%:%.2f%%)\n", s->n_sgltn[0], s->n_sgltn[1], (float)s->n_sgltn[0] / s->n_pair_all[0] * 100.0, (float)s->n_sgltn[1] / s->n_pair_all[1] * 100.0);
	printf("%lld + %lld with mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
	printf("%lld + %lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
	free(s);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}
