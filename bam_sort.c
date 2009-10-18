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
	uint64_t pos, idx;
	bam1_t *b;
} heap1_t;

#define __pos_cmp(a, b) ((a).pos > (b).pos || ((a).pos == (b).pos && ((a).i > (b).i || ((a).i == (b).i && (a).idx > (b).idx))))

static inline int heap_lt(const heap1_t a, const heap1_t b)
{
	if (g_is_by_qname) {
		int t;
		if (a.b == 0 || b.b == 0) return a.b == 0? 1 : 0;
		t = strnum_cmp(bam1_qname(a.b), bam1_qname(b.b));
		return (t > 0 || (t == 0 && __pos_cmp(a, b)));
	} else return __pos_cmp(a, b);
}

KSORT_INIT(heap, heap1_t, heap_lt)

static void swap_header_text(bam_header_t *h1, bam_header_t *h2)
{
	int tempi;
	char *temps;
	tempi = h1->l_text, h1->l_text = h2->l_text, h2->l_text = tempi;
	temps = h1->text, h1->text = h2->text, h2->text = temps;
}

/*!
  @abstract    Merge multiple sorted BAM.
  @param  is_by_qname whether to sort by query name
  @param  out  output BAM file name
  @param  headers  name of SAM file from which to copy '@' header lines,
                   or NULL to copy them from the first file to be merged
  @param  n    number of files to be merged
  @param  fn   names of files to be merged

  @discussion Padding information may NOT correctly maintained. This
  function is NOT thread safe.
 */
void bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int add_RG)
{
	bamFile fpout, *fp;
	heap1_t *heap;
	bam_header_t *hout = 0;
	bam_header_t *hheaders = NULL;
	int i, j, *RG_len = 0;
	uint64_t idx = 0;
	char **RG = 0;

	if (headers) {
		tamFile fpheaders = sam_open(headers);
		if (fpheaders == 0) {
			fprintf(stderr, "[bam_merge_core] Cannot open file `%s'. Continue anyway.\n", headers);
		} else {
			hheaders = sam_header_read(fpheaders);
			sam_close(fpheaders);
		}
	}

	g_is_by_qname = by_qname;
	fp = (bamFile*)calloc(n, sizeof(bamFile));
	heap = (heap1_t*)calloc(n, sizeof(heap1_t));
	// prepare RG tag
	if (add_RG) {
		RG = (char**)calloc(n, sizeof(void*));
		RG_len = (int*)calloc(n, sizeof(int));
		for (i = 0; i != n; ++i) {
			int l = strlen(fn[i]);
			const char *s = fn[i];
			if (l > 4 && strcmp(s + l - 4, ".bam") == 0) l -= 4;
			for (j = l - 1; j >= 0; --j) if (s[j] == '/') break;
			++j; l -= j;
			RG[i] = calloc(l + 1, 1);
			RG_len[i] = l;
			strncpy(RG[i], s + j, l);
		}
	}
	// read the first
	for (i = 0; i != n; ++i) {
		heap1_t *h;
		bam_header_t *hin;
		fp[i] = bam_open(fn[i], "r");
		if (fp[i] == 0) {
			int j;
			fprintf(stderr, "[bam_merge_core] fail to open file %s\n", fn[i]);
			for (j = 0; j < i; ++j) bam_close(fp[j]);
			free(fp); free(heap);
			// FIXME: possible memory leak
			return;
		}
		hin = bam_header_read(fp[i]);
		if (i == 0) { // the first SAM
			hout = hin;
			if (hheaders) {
				// If the text headers to be swapped in include any @SQ headers,
				// check that they are consistent with the existing binary list
				// of reference information.
				if (hheaders->n_targets > 0) {
					if (hout->n_targets != hheaders->n_targets)
						fprintf(stderr, "[bam_merge_core] number of @SQ headers in `%s' differs from number of target sequences", headers);
					for (j = 0; j < hout->n_targets; ++j)
						if (strcmp(hout->target_name[j], hheaders->target_name[j]) != 0)
							fprintf(stderr, "[bam_merge_core] @SQ header '%s' in '%s' differs from target sequence", hheaders->target_name[j], headers);
				}
				swap_header_text(hout, hheaders);
				bam_header_destroy(hheaders);
				hheaders = NULL;
			}
		} else { // validate multiple baf
			if (hout->n_targets != hin->n_targets) {
				fprintf(stderr, "[bam_merge_core] file '%s' has different number of target sequences. Abort!\n", fn[i]);
				exit(1);
			}
			for (j = 0; j < hout->n_targets; ++j) {
				if (strcmp(hout->target_name[j], hin->target_name[j])) {
					fprintf(stderr, "[bam_merge_core] different target sequence name: '%s' != '%s' in file '%s'. Abort!\n",
							hout->target_name[j], hin->target_name[j], fn[i]);
					exit(1);
				}
			}
			bam_header_destroy(hin);
		}
		h = heap + i;
		h->i = i;
		h->b = (bam1_t*)calloc(1, sizeof(bam1_t));
		if (bam_read1(fp[i], h->b) >= 0) {
			h->pos = ((uint64_t)h->b->core.tid<<32) | (uint32_t)h->b->core.pos<<1 | bam1_strand(h->b);
			h->idx = idx++;
		}
		else h->pos = HEAP_EMPTY;
	}
	fpout = strcmp(out, "-")? bam_open(out, "w") : bam_dopen(fileno(stdout), "w");
	assert(fpout);
	bam_header_write(fpout, hout);
	bam_header_destroy(hout);

	ks_heapmake(heap, n, heap);
	while (heap->pos != HEAP_EMPTY) {
		bam1_t *b = heap->b;
		if (add_RG && bam_aux_get(b, "RG") == 0)
			bam_aux_append(b, "RG", 'Z', RG_len[heap->i] + 1, (uint8_t*)RG[heap->i]);
		bam_write1_core(fpout, &b->core, b->data_len, b->data);
		if ((j = bam_read1(fp[heap->i], b)) >= 0) {
			heap->pos = ((uint64_t)b->core.tid<<32) | (uint32_t)b->core.pos<<1 | bam1_strand(b);
			heap->idx = idx++;
		} else if (j == -1) {
			heap->pos = HEAP_EMPTY;
			free(heap->b->data); free(heap->b);
			heap->b = 0;
		} else fprintf(stderr, "[bam_merge_core] '%s' is truncated. Continue anyway.\n", fn[heap->i]);
		ks_heapadjust(heap, 0, n, heap);
	}

	if (add_RG) {
		for (i = 0; i != n; ++i) free(RG[i]);
		free(RG); free(RG_len);
	}
	for (i = 0; i != n; ++i) bam_close(fp[i]);
	bam_close(fpout);
	free(fp); free(heap);
}
int bam_merge(int argc, char *argv[])
{
	int c, is_by_qname = 0, add_RG = 0;
	char *fn_headers = NULL;

	while ((c = getopt(argc, argv, "h:nr")) >= 0) {
		switch (c) {
		case 'r': add_RG = 1; break;
		case 'h': fn_headers = strdup(optarg); break;
		case 'n': is_by_qname = 1; break;
		}
	}
	if (optind + 2 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools merge [-nr] [-h inh.sam] <out.bam> <in1.bam> <in2.bam> [...]\n\n");
		fprintf(stderr, "Options: -n       sort by read names\n");
		fprintf(stderr, "         -r       attach RG tag (inferred from file names)\n");
		fprintf(stderr, "         -h FILE  copy the header in FILE to <out.bam> [in1.bam]\n\n");
		fprintf(stderr, "Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users\n");
		fprintf(stderr, "      must provide the correct header with -h, or uses Picard which properly maintains\n");
		fprintf(stderr, "      the header dictionary in merging.\n\n");
		return 1;
	}
	bam_merge_core(is_by_qname, argv[optind], fn_headers, argc - optind - 1, argv + optind + 1, add_RG);
	free(fn_headers);
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

static void sort_blocks(int n, int k, bam1_p *buf, const char *prefix, const bam_header_t *h, int is_stdout)
{
	char *name;
	int i;
	bamFile fp;
	ks_mergesort(sort, k, buf, 0);
	name = (char*)calloc(strlen(prefix) + 20, 1);
	if (n >= 0) sprintf(name, "%s.%.4d.bam", prefix, n);
	else sprintf(name, "%s.bam", prefix);
	fp = is_stdout? bam_dopen(fileno(stdout), "w") : bam_open(name, "w");
	if (fp == 0) {
		fprintf(stderr, "[sort_blocks] fail to create file %s.\n", name);
		free(name);
		// FIXME: possible memory leak
		return;
	}
	free(name);
	bam_header_write(fp, h);
	for (i = 0; i < k; ++i)
		bam_write1_core(fp, &buf[i]->core, buf[i]->data_len, buf[i]->data);
	bam_close(fp);
}

/*!
  @abstract Sort an unsorted BAM file based on the chromosome order
  and the leftmost position of an alignment

  @param  is_by_qname whether to sort by query name
  @param  fn       name of the file to be sorted
  @param  prefix   prefix of the output and the temporary files; upon
	                   sucessess, prefix.bam will be written.
  @param  max_mem  approxiate maximum memory (very inaccurate)

  @discussion It may create multiple temporary subalignment files
  and then merge them by calling bam_merge_core(). This function is
  NOT thread safe.
 */
void bam_sort_core_ext(int is_by_qname, const char *fn, const char *prefix, size_t max_mem, int is_stdout)
{
	int n, ret, k, i;
	size_t mem;
	bam_header_t *header;
	bamFile fp;
	bam1_t *b, **buf;

	g_is_by_qname = is_by_qname;
	n = k = 0; mem = 0;
	fp = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[bam_sort_core] fail to open file %s\n", fn);
		return;
	}
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
			sort_blocks(n++, k, buf, prefix, header, is_stdout);
			mem = 0; k = 0;
		}
	}
	if (ret != -1)
		fprintf(stderr, "[bam_sort_core] truncated file. Continue anyway.\n");
	if (n == 0) sort_blocks(-1, k, buf, prefix, header, is_stdout);
	else { // then merge
		char **fns, *fnout;
		fprintf(stderr, "[bam_sort_core] merging from %d files...\n", n+1);
		sort_blocks(n++, k, buf, prefix, header, is_stdout);
		fnout = (char*)calloc(strlen(prefix) + 20, 1);
		if (is_stdout) sprintf(fnout, "-");
		else sprintf(fnout, "%s.bam", prefix);
		fns = (char**)calloc(n, sizeof(char*));
		for (i = 0; i < n; ++i) {
			fns[i] = (char*)calloc(strlen(prefix) + 20, 1);
			sprintf(fns[i], "%s.%.4d.bam", prefix, i);
		}
		bam_merge_core(is_by_qname, fnout, 0, n, fns, 0);
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

void bam_sort_core(int is_by_qname, const char *fn, const char *prefix, size_t max_mem)
{
	bam_sort_core_ext(is_by_qname, fn, prefix, max_mem, 0);
}

int bam_sort(int argc, char *argv[])
{
	size_t max_mem = 500000000;
	int c, is_by_qname = 0, is_stdout = 0;
	while ((c = getopt(argc, argv, "nom:")) >= 0) {
		switch (c) {
		case 'o': is_stdout = 1; break;
		case 'n': is_by_qname = 1; break;
		case 'm': max_mem = atol(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: samtools sort [-on] [-m <maxMem>] <in.bam> <out.prefix>\n");
		return 1;
	}
	bam_sort_core_ext(is_by_qname, argv[optind], argv[optind+1], max_mem, is_stdout);
	return 0;
}
