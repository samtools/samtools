#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "kstring.h"
#include "bgzf.h"
#include "bam.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

typedef struct {
	bamFile fp;
	bam_iter_t iter;
	int min_mapQ;
} aux_t;

static int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data;
	int ret = bam_iter_read(aux->fp, aux->iter, b);
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	return ret;
}

int main_bedcov(int argc, char *argv[])
{
	extern void bam_init_header_hash(bam_header_t*);
	gzFile fp;
	kstring_t str;
	kstream_t *ks;
	bam_index_t **idx;
	bam_header_t *h = 0;
	aux_t **aux;
	int *n_plp, dret, i, n, c, min_mapQ = 0;
	int64_t *cnt;
	const bam_pileup1_t **plp;

	while ((c = getopt(argc, argv, "Q:")) >= 0) {
		switch (c) {
		case 'Q': min_mapQ = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: samtools bedcov <in.bed> <in1.bam> [...]\n");
		return 1;
	}
	memset(&str, 0, sizeof(kstring_t));
	n = argc - optind - 1;
	aux = calloc(n, sizeof(void*));
	idx = calloc(n, sizeof(void*));
	for (i = 0; i < n; ++i) {
		aux[i] = calloc(1, sizeof(aux_t));
		aux[i]->min_mapQ = min_mapQ;
		aux[i]->fp = bam_open(argv[i+optind+1], "r");
		idx[i] = bam_index_load(argv[i+optind+1]);
		if (aux[i]->fp == 0 || idx[i] == 0) {
			fprintf(stderr, "ERROR: fail to open index BAM file '%s'\n", argv[i+optind+1]);
			return 2;
		}
		bgzf_set_cache_size(aux[i]->fp, 20);
		if (i == 0) h = bam_header_read(aux[0]->fp);
	}
	bam_init_header_hash(h);
	cnt = calloc(n, 8);

	fp = gzopen(argv[optind], "rb");
	ks = ks_init(fp);
	n_plp = calloc(n, sizeof(int));
	plp = calloc(n, sizeof(void*));
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *p, *q;
		int tid, beg, end, pos;
		bam_mplp_t mplp;

		for (p = q = str.s; *p && *p != '\t'; ++p);
		if (*p != '\t') goto bed_error;
		*p = 0; tid = bam_get_tid(h, q); *p = '\t';
		if (tid < 0) goto bed_error;
		for (q = p = p + 1; isdigit(*p); ++p);
		if (*p != '\t') goto bed_error;
		*p = 0; beg = atoi(q); *p = '\t';
		for (q = p = p + 1; isdigit(*p); ++p);
		if (*p == '\t' || *p == 0) {
			int c = *p;
			*p = 0; end = atoi(q); *p = c;
		} else goto bed_error;

		for (i = 0; i < n; ++i) {
			if (aux[i]->iter) bam_iter_destroy(aux[i]->iter);
			aux[i]->iter = bam_iter_query(idx[i], tid, beg, end);
		}
		mplp = bam_mplp_init(n, read_bam, (void**)aux);
		bam_mplp_set_maxcnt(mplp, 64000);
		memset(cnt, 0, 8 * n);
		while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0)
			if (pos >= beg && pos < end)
				for (i = 0; i < n; ++i) cnt[i] += n_plp[i];
		for (i = 0; i < n; ++i) {
			kputc('\t', &str);
			kputl(cnt[i], &str);
		}
		puts(str.s);
		bam_mplp_destroy(mplp);
		continue;

bed_error:
		fprintf(stderr, "Errors in BED line '%s'\n", str.s);
	}
	free(n_plp); free(plp);
	ks_destroy(ks);
	gzclose(fp);

	free(cnt);
	for (i = 0; i < n; ++i) {
		if (aux[i]->iter) bam_iter_destroy(aux[i]->iter);
		bam_index_destroy(idx[i]);
		bam_close(aux[i]->fp);
		free(aux[i]);
	}
	bam_header_destroy(h);
	free(aux); free(idx);
	free(str.s);
	return 0;
}
