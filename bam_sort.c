#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bam.h"
#include "ksort.h"

static int g_is_by_qname = 0;

static inline int strnum_cmp(const char *a, const char *b)
{
	char *pa, *pb;
	pa = (char*)a; pb = (char*)b;
	while (*pa && *pb) {
		if (isdigit(*pa) && isdigit(*pb)) {
			long ai, bi;
			ai = strtol(pa, &pa, 10);
			bi = strtol(pb, &pb, 10);
			if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
		} else {
			if (*pa != *pb) break;
			++pa; ++pb;
		}
	}
	if (*pa == *pb)
		return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
	return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

#define HEAP_EMPTY 0xffffffffffffffffull

typedef struct {
	int i;
	uint64_t pos;
	bam1_t *b;
} heap1_t;

static inline int heap_lt(const heap1_t a, const heap1_t b)
{
	if (g_is_by_qname) {
		int t = strnum_cmp(bam1_qname(a.b), bam1_qname(b.b));
		return (t > 0 || (t == 0 && a.pos > b.pos));
	} else return (a.pos > b.pos);
}

KSORT_INIT(heap, heap1_t, heap_lt)

void bam_merge_core(int by_qname, const char *out, int n, char * const *fn)
{
	bamFile fpout, *fp;
	heap1_t *heap;
	bam_header_t *hout = 0;
	int i, j;

	g_is_by_qname = by_qname;
	fp = (bamFile*)calloc(n, sizeof(bamFile));
	heap = (heap1_t*)calloc(n, sizeof(heap1_t));
	for (i = 0; i != n; ++i) {
		heap1_t *h;
		bam_header_t *hin;
		assert(fp[i] = bam_open(fn[i], "r"));
		hin = bam_header_read(fp[i]);
		if (i == 0) hout = hin;
		else { // validate multiple baf
			if (hout->n_targets != hin->n_targets) {
				fprintf(stderr, "[bam_merge_core] file '%s' has different number of target sequences. Abort!\n", fn[i]);
				abort();
			}
			for (j = 0; j < hout->n_targets; ++j) {
				if (strcmp(hout->target_name[j], hin->target_name[j]) || hout->target_len[j] != hin->target_len[j]) {
					fprintf(stderr, "[bam_merge_core] file '%s' has a different target sequence. Abort!\n", fn[i]);
					abort();
				}
			}
			bam_header_destroy(hin);
		}
		h = heap + i;
		h->i = i;
		h->b = (bam1_t*)calloc(1, sizeof(bam1_t));
		if (bam_read1(fp[i], h->b) >= 0)
			h->pos = ((uint64_t)h->b->core.tid<<32) | (uint32_t)h->b->core.pos<<1 | bam1_strand(h->b);
		else h->pos = HEAP_EMPTY;
	}
	fpout = strcmp(out, "-")? bam_open(out, "w") : bam_dopen(fileno(stdout), "w");
	assert(fpout);
	bam_header_write(fpout, hout);
	bam_header_destroy(hout);

	ks_heapmake(heap, n, heap);
	while (heap->pos != HEAP_EMPTY) {
		bam1_t *b = heap->b;
		bam_write1_core(fpout, &b->core, b->data_len, b->data);
		if ((j = bam_read1(fp[heap->i], b)) >= 0)
			heap->pos = ((uint64_t)b->core.tid<<32) | (uint32_t)b->core.pos<<1 | bam1_strand(b);
		else if (j == -1) heap->pos = HEAP_EMPTY;
		else fprintf(stderr, "[bam_merge_core] '%s' is truncated. Continue anyway.\n", fn[heap->i]);
		ks_heapadjust(heap, 0, n, heap);
	}

	for (i = 0; i != n; ++i) {
		bam_close(fp[i]);
		free(heap[i].b->data);
		free(heap[i].b);
	}
	bam_close(fpout);
	free(fp); free(heap);
}
int bam_merge(int argc, char *argv[])
{
	int c, is_by_qname = 0;
	while ((c = getopt(argc, argv, "n")) >= 0) {
		switch (c) {
		case 'n': is_by_qname = 1; break;
		}
	}
	if (optind + 3 >= argc) {
		fprintf(stderr, "Usage: samtools merge [-n] <out.bam> <in1.bam> <in2.bam> [...]\n");
		return 1;
	}
	bam_merge_core(is_by_qname, argv[optind], argc - optind - 1, argv + optind + 1);
	return 0;
}

typedef bam1_t *bam1_p;

static inline int bam1_lt(const bam1_p a, const bam1_p b)
{
	if (g_is_by_qname) {
		int t = strnum_cmp(bam1_qname(a), bam1_qname(b));
		return (t < 0 || (t == 0 && (((uint64_t)a->core.tid<<32|a->core.pos) < ((uint64_t)b->core.tid<<32|b->core.pos))));
	} else return (((uint64_t)a->core.tid<<32|a->core.pos) < ((uint64_t)b->core.tid<<32|b->core.pos));
}
KSORT_INIT(sort, bam1_p, bam1_lt)

static void sort_blocks(int n, int k, bam1_p *buf, const char *prefix, const bam_header_t *h)
{
	char *name;
	int i;
	bamFile fp;
	ks_mergesort(sort, k, buf, 0);
	name = (char*)calloc(strlen(prefix) + 20, 1);
	if (n >= 0) sprintf(name, "%s.%.4d.bam", prefix, n);
	else sprintf(name, "%s.bam", prefix);
	assert(fp = bam_open(name, "w"));
	free(name);
	bam_header_write(fp, h);
	for (i = 0; i < k; ++i)
		bam_write1_core(fp, &buf[i]->core, buf[i]->data_len, buf[i]->data);
	bam_close(fp);
}

void bam_sort_core(int is_by_qname, const char *fn, const char *prefix, size_t max_mem)
{
	int n, ret, k, i;
	size_t mem;
	bam_header_t *header;
	bamFile fp;
	bam1_t *b, **buf;

	g_is_by_qname = is_by_qname;
	n = k = 0; mem = 0;
	fp = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);
	buf = (bam1_t**)calloc(max_mem / BAM_CORE_SIZE, sizeof(bam1_t*));
	// write sub files
	for (;;) {
		if (buf[k] == 0) buf[k] = (bam1_t*)calloc(1, sizeof(bam1_t));
		b = buf[k];
		if ((ret = bam_read1(fp, b)) < 0) break;
		mem += ret;
		++k;
		if (mem >= max_mem) {
			sort_blocks(n++, k, buf, prefix, header);
			mem = 0; k = 0;
		}
	}
	if (ret != -1)
		fprintf(stderr, "[bam_sort_core] truncated file. Continue anyway.\n");
	if (n == 0) sort_blocks(-1, k, buf, prefix, header);
	else { // then merge
		char **fns, *fnout;
		fprintf(stderr, "[bam_sort_core] merging from %d files...\n", n+1);
		sort_blocks(n++, k, buf, prefix, header);
		fnout = (char*)calloc(strlen(prefix) + 20, 1);
		sprintf(fnout, "%s.bam", prefix);
		fns = (char**)calloc(n, sizeof(char*));
		for (i = 0; i < n; ++i) {
			fns[i] = (char*)calloc(strlen(prefix) + 20, 1);
			sprintf(fns[i], "%s.%.4d.bam", prefix, i);
		}
		bam_merge_core(is_by_qname, fnout, n, fns);
		free(fnout);
		for (i = 0; i < n; ++i) {
			unlink(fns[i]);
			free(fns[i]);
		}
		free(fns);
	}
	for (k = 0; k < max_mem / BAM_CORE_SIZE; ++k) {
		if (buf[k]) {
			free(buf[k]->data);
			free(buf[k]);
		}
	}
	free(buf);
	bam_header_destroy(header);
	bam_close(fp);
}

int bam_sort(int argc, char *argv[])
{
	size_t max_mem = 500000000;
	int c, is_by_qname = 0;
	while ((c = getopt(argc, argv, "nm:")) >= 0) {
		switch (c) {
		case 'n': is_by_qname = 1; break;
		case 'm': max_mem = atol(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: samtools sort [-n] [-m <maxMem>] <in.baf> <out.prefix>\n");
		return 1;
	}
	bam_sort_core(is_by_qname, argv[optind], argv[optind+1], max_mem);
	return 0;
}
