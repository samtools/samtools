#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "bam.h"

typedef struct {
	bamFile fp;
	bam_iter_t iter;
} aux_t;

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

static int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data;
	return aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
}

#ifdef BAM2BED_MAIN
int main(int argc, char *argv[])
#else
int main_depth(int argc, char *argv[])
#endif
{
	int i, n, tid, beg, end, pos, *n_plp;
	const bam_pileup1_t **plp;
	char *reg = 0;
	void *bed = 0;
	bam_header_t *h = 0;
	aux_t **data;
	bam_mplp_t mplp;

	while ((n = getopt(argc, argv, "r:b:")) >= 0) {
		switch (n) {
			case 'r': reg = strdup(optarg); break;
			case 'b': bed = bed_read(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: bam2depth [-r reg] [-b in.bed] <in1.bam> [...]\n");
		return 1;
	}
	n = argc - optind;
	data = calloc(n, sizeof(void*));
	tid = -1; beg = 0; end = 1<<30;
	for (i = 0; i < n; ++i) {
		bam_header_t *htmp;
		data[i] = calloc(1, sizeof(aux_t));
		data[i]->fp = bam_open(argv[optind+i], "r");
		htmp = bam_header_read(data[i]->fp);
		if (i == 0) {
			h = htmp;
			if (reg) bam_parse_region(h, reg, &tid, &beg, &end);
		} else bam_header_destroy(htmp);
		if (tid >= 0) {
			bam_index_t *idx = bam_index_load(argv[optind+i]);
			data[i]->iter = bam_iter_query(idx, tid, beg, end);
			bam_index_destroy(idx);
		}
	}
	mplp = bam_mplp_init(n, read_bam, (void**)data);
	n_plp = calloc(n, sizeof(int));
	plp = calloc(n, sizeof(void*));
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
		if (pos < beg || pos >= end) continue;
		if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue;
		printf("%s\t%d", h->target_name[tid], pos+1);
		for (i = 0; i < n; ++i) {
			int j, m = 0;
			const bam_pileup1_t *p = plp[i];
			for (j = 0; j < n_plp[i]; ++j)
				if (p->is_del || p->is_refskip) ++m;
			printf("\t%d", n_plp[i] - m);
		}
		putchar('\n');
	}
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);
	bam_header_destroy(h);
	for (i = 0; i < n; ++i) {
		bam_close(data[i]->fp);
		if (data[i]->iter) bam_iter_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data);
	if (reg) free(reg);
	if (bed) bed_destroy(bed);
	return 0;
}
