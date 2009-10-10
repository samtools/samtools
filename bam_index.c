#include <ctype.h>
#include <assert.h>
#include "bam.h"
#include "khash.h"
#include "ksort.h"
#include "bam_endian.h"
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

/*!
  @header

  Alignment indexing. Before indexing, BAM must be sorted based on the
  leftmost coordinate of alignments. In indexing, BAM uses two indices:
  a UCSC binning index and a simple linear index. The binning index is
  efficient for alignments spanning long distance, while the auxiliary
  linear index helps to reduce unnecessary seek calls especially for
  short alignments.

  The UCSC binning scheme was suggested by Richard Durbin and Lincoln
  Stein and is explained by Kent et al. (2002). In this scheme, each bin
  represents a contiguous genomic region which can be fully contained in
  another bin; each alignment is associated with a bin which represents
  the smallest region containing the entire alignment. The binning
  scheme is essentially another representation of R-tree. A distinct bin
  uniquely corresponds to a distinct internal node in a R-tree. Bin A is
  a child of Bin B if region A is contained in B.

  In BAM, each bin may span 2^29, 2^26, 2^23, 2^20, 2^17 or 2^14 bp. Bin
  0 spans a 512Mbp region, bins 1-8 span 64Mbp, 9-72 8Mbp, 73-584 1Mbp,
  585-4680 128Kbp and bins 4681-37449 span 16Kbp regions. If we want to
  find the alignments overlapped with a region [rbeg,rend), we need to
  calculate the list of bins that may be overlapped the region and test
  the alignments in the bins to confirm the overlaps. If the specified
  region is short, typically only a few alignments in six bins need to
  be retrieved. The overlapping alignments can be quickly fetched.

 */

#define BAM_MIN_CHUNK_GAP 32768
// 1<<14 is the size of minimum bin.
#define BAM_LIDX_SHIFT    14

typedef struct {
	uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(off, pair64_t, pair64_lt)

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} bam_binlist_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} bam_lidx_t;

KHASH_MAP_INIT_INT(i, bam_binlist_t)

struct __bam_index_t {
	int32_t n;
	khash_t(i) **index;
	bam_lidx_t *index2;
};

// requirement: len <= LEN_MASK
static inline void insert_offset(khash_t(i) *h, int bin, uint64_t beg, uint64_t end)
{
	khint_t k;
	bam_binlist_t *l;
	int ret;
	k = kh_put(i, h, bin, &ret);
	l = &kh_value(h, k);
	if (ret) { // not present
		l->m = 1; l->n = 0;
		l->list = (pair64_t*)calloc(l->m, 16);
	}
	if (l->n == l->m) {
		l->m <<= 1;
		l->list = (pair64_t*)realloc(l->list, l->m * 16);
	}
	l->list[l->n].u = beg; l->list[l->n++].v = end;
}

static inline void insert_offset2(bam_lidx_t *index2, bam1_t *b, uint64_t offset)
{
	int i, beg, end;
	beg = b->core.pos >> BAM_LIDX_SHIFT;
	end = (bam_calend(&b->core, bam1_cigar(b)) - 1) >> BAM_LIDX_SHIFT;
	if (index2->m < end + 1) {
		int old_m = index2->m;
		index2->m = end + 1;
		kroundup32(index2->m);
		index2->offset = (uint64_t*)realloc(index2->offset, index2->m * 8);
		memset(index2->offset + old_m, 0, 8 * (index2->m - old_m));
	}
	for (i = beg + 1; i <= end; ++i)
		if (index2->offset[i] == 0) index2->offset[i] = offset;
	index2->n = end + 1;
}

static void merge_chunks(bam_index_t *idx)
{
#if defined(BAM_TRUE_OFFSET) || defined(BAM_VIRTUAL_OFFSET16)
	khash_t(i) *index;
	int i, l, m;
	khint_t k;
	for (i = 0; i < idx->n; ++i) {
		index = idx->index[i];
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			bam_binlist_t *p;
			if (!kh_exist(index, k)) continue;
			p = &kh_value(index, k);
			m = 0;
			for (l = 1; l < p->n; ++l) {
#ifdef BAM_TRUE_OFFSET
				if (p->list[m].v + BAM_MIN_CHUNK_GAP > p->list[l].u) p->list[m].v = p->list[l].v;
#else
				if (p->list[m].v>>16 == p->list[l].u>>16) p->list[m].v = p->list[l].v;
#endif
				else p->list[++m] = p->list[l];
			} // ~for(l)
			p->n = m + 1;
		} // ~for(k)
	} // ~for(i)
#endif // defined(BAM_TRUE_OFFSET) || defined(BAM_BGZF)
}

bam_index_t *bam_index_core(bamFile fp)
{
	bam1_t *b;
	bam_header_t *h;
	int i, ret;
	bam_index_t *idx;
	uint32_t last_bin, save_bin;
	int32_t last_coor, last_tid, save_tid;
	bam1_core_t *c;
	uint64_t save_off, last_off;

	idx = (bam_index_t*)calloc(1, sizeof(bam_index_t));
	b = (bam1_t*)calloc(1, sizeof(bam1_t));
	h = bam_header_read(fp);
	c = &b->core;

	idx->n = h->n_targets;
	bam_header_destroy(h);
	idx->index = (khash_t(i)**)calloc(idx->n, sizeof(void*));
	for (i = 0; i < idx->n; ++i) idx->index[i] = kh_init(i);
	idx->index2 = (bam_lidx_t*)calloc(idx->n, sizeof(bam_lidx_t));

	save_bin = save_tid = last_tid = last_bin = 0xffffffffu;
	save_off = last_off = bam_tell(fp); last_coor = 0xffffffffu;
	while ((ret = bam_read1(fp, b)) >= 0) {
		if (last_tid != c->tid) { // change of chromosomes
			last_tid = c->tid;
			last_bin = 0xffffffffu;
		} else if (last_coor > c->pos) {
			fprintf(stderr, "[bam_index_core] the alignment is not sorted (%s): %u > %u in %d-th chr\n",
					bam1_qname(b), last_coor, c->pos, c->tid+1);
			exit(1);
		}
		if (b->core.tid >= 0 && b->core.bin < 4681) insert_offset2(&idx->index2[b->core.tid], b, last_off);
		if (c->bin != last_bin) { // then possibly write the binning index
			if (save_bin != 0xffffffffu) // save_bin==0xffffffffu only happens to the first record
				insert_offset(idx->index[save_tid], save_bin, save_off, last_off);
			save_off = last_off;
			save_bin = last_bin = c->bin;
			save_tid = c->tid;
			if (save_tid < 0) break;
		}
		if (bam_tell(fp) <= last_off) {
			fprintf(stderr, "[bam_index_core] bug in BGZF/RAZF: %llx < %llx\n",
					(unsigned long long)bam_tell(fp), (unsigned long long)last_off);
			exit(1);
		}
		last_off = bam_tell(fp);
		last_coor = b->core.pos;
	}
	if (save_tid >= 0) insert_offset(idx->index[save_tid], save_bin, save_off, bam_tell(fp));
	merge_chunks(idx);
	if (ret < -1) fprintf(stderr, "[bam_index_core] truncated file? Continue anyway. (%d)\n", ret);
	free(b->data); free(b);
	return idx;
}

void bam_index_destroy(bam_index_t *idx)
{
	khint_t k;
	int i;
	if (idx == 0) return;
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index = idx->index[i];
		bam_lidx_t *index2 = idx->index2 + i;
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			if (kh_exist(index, k))
				free(kh_value(index, k).list);
		}
		kh_destroy(i, index);
		free(index2->offset);
	}
	free(idx->index); free(idx->index2);
	free(idx);
}

void bam_index_save(const bam_index_t *idx, FILE *fp)
{
	int32_t i, size;
	khint_t k;
	fwrite("BAI\1", 1, 4, fp);
	if (bam_is_be) {
		uint32_t x = idx->n;
		fwrite(bam_swap_endian_4p(&x), 4, 1, fp);
	} else fwrite(&idx->n, 4, 1, fp);
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index = idx->index[i];
		bam_lidx_t *index2 = idx->index2 + i;
		// write binning index
		size = kh_size(index);
		if (bam_is_be) { // big endian
			uint32_t x = size;
			fwrite(bam_swap_endian_4p(&x), 4, 1, fp);
		} else fwrite(&size, 4, 1, fp);
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			if (kh_exist(index, k)) {
				bam_binlist_t *p = &kh_value(index, k);
				if (bam_is_be) { // big endian
					uint32_t x;
					x = kh_key(index, k); fwrite(bam_swap_endian_4p(&x), 4, 1, fp);
					x = p->n; fwrite(bam_swap_endian_4p(&x), 4, 1, fp);
					for (x = 0; (int)x < p->n; ++x) {
						bam_swap_endian_8p(&p->list[x].u);
						bam_swap_endian_8p(&p->list[x].v);
					}
					fwrite(p->list, 16, p->n, fp);
					for (x = 0; (int)x < p->n; ++x) {
						bam_swap_endian_8p(&p->list[x].u);
						bam_swap_endian_8p(&p->list[x].v);
					}
				} else {
					fwrite(&kh_key(index, k), 4, 1, fp);
					fwrite(&p->n, 4, 1, fp);
					fwrite(p->list, 16, p->n, fp);
				}
			}
		}
		// write linear index (index2)
		if (bam_is_be) {
			int x = index2->n;
			fwrite(bam_swap_endian_4p(&x), 4, 1, fp);
		} else fwrite(&index2->n, 4, 1, fp);
		if (bam_is_be) { // big endian
			int x;
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
			fwrite(index2->offset, 8, index2->n, fp);
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
		} else fwrite(index2->offset, 8, index2->n, fp);
	}
	fflush(fp);
}

static bam_index_t *bam_index_load_core(FILE *fp)
{
	int i;
	char magic[4];
	bam_index_t *idx;
	if (fp == 0) {
		fprintf(stderr, "[bam_index_load_core] fail to load index.\n");
		return 0;
	}
	fread(magic, 1, 4, fp);
	if (strncmp(magic, "BAI\1", 4)) {
		fprintf(stderr, "[bam_index_load] wrong magic number.\n");
		fclose(fp);
		return 0;
	}
	idx = (bam_index_t*)calloc(1, sizeof(bam_index_t));	
	fread(&idx->n, 4, 1, fp);
	if (bam_is_be) bam_swap_endian_4p(&idx->n);
	idx->index = (khash_t(i)**)calloc(idx->n, sizeof(void*));
	idx->index2 = (bam_lidx_t*)calloc(idx->n, sizeof(bam_lidx_t));
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index;
		bam_lidx_t *index2 = idx->index2 + i;
		uint32_t key, size;
		khint_t k;
		int j, ret;
		bam_binlist_t *p;
		index = idx->index[i] = kh_init(i);
		// load binning index
		fread(&size, 4, 1, fp);
		if (bam_is_be) bam_swap_endian_4p(&size);
		for (j = 0; j < (int)size; ++j) {
			fread(&key, 4, 1, fp);
			if (bam_is_be) bam_swap_endian_4p(&key);
			k = kh_put(i, index, key, &ret);
			p = &kh_value(index, k);
			fread(&p->n, 4, 1, fp);
			if (bam_is_be) bam_swap_endian_4p(&p->n);
			p->m = p->n;
			p->list = (pair64_t*)malloc(p->m * 16);
			fread(p->list, 16, p->n, fp);
			if (bam_is_be) {
				int x;
				for (x = 0; x < p->n; ++x) {
					bam_swap_endian_8p(&p->list[x].u);
					bam_swap_endian_8p(&p->list[x].v);
				}
			}
		}
		// load linear index
		fread(&index2->n, 4, 1, fp);
		if (bam_is_be) bam_swap_endian_4p(&index2->n);
		index2->m = index2->n;
		index2->offset = (uint64_t*)calloc(index2->m, 8);
		fread(index2->offset, index2->n, 8, fp);
		if (bam_is_be)
			for (j = 0; j < index2->n; ++j) bam_swap_endian_8p(&index2->offset[j]);
	}
	return idx;
}

bam_index_t *bam_index_load_local(const char *_fn)
{
	FILE *fp;
	char *fnidx, *fn;

	if (strstr(_fn, "ftp://") == _fn || strstr(_fn, "http://") == _fn) {
		const char *p;
		int l = strlen(_fn);
		for (p = _fn + l - 1; p >= _fn; --p)
			if (*p == '/') break;
		fn = strdup(p + 1);
	} else fn = strdup(_fn);
	fnidx = (char*)calloc(strlen(fn) + 5, 1);
	strcpy(fnidx, fn); strcat(fnidx, ".bai");
	fp = fopen(fnidx, "r");
	if (fp == 0) { // try "{base}.bai"
		char *s = strstr(fn, "bam");
		if (s == fn + strlen(fn) - 3) {
			strcpy(fnidx, fn);
			fnidx[strlen(fn)-1] = 'i';
			fp = fopen(fnidx, "r");
		}
	}
	free(fnidx); free(fn);
	if (fp) {
		bam_index_t *idx = bam_index_load_core(fp);
		fclose(fp);
		return idx;
	} else return 0;
}

#ifdef _USE_KNETFILE
static void download_from_remote(const char *url)
{
	const int buf_size = 1 * 1024 * 1024;
	char *fn;
	FILE *fp;
	uint8_t *buf;
	knetFile *fp_remote;
	int l;
	if (strstr(url, "ftp://") != url && strstr(url, "http://") != url) return;
	l = strlen(url);
	for (fn = (char*)url + l - 1; fn >= url; --fn)
		if (*fn == '/') break;
	++fn; // fn now points to the file name
	fp_remote = knet_open(url, "r");
	if (fp_remote == 0) {
		fprintf(stderr, "[download_from_remote] fail to open remote file.\n");
		return;
	}
	if ((fp = fopen(fn, "w")) == 0) {
		fprintf(stderr, "[download_from_remote] fail to create file in the working directory.\n");
		knet_close(fp_remote);
		return;
	}
	buf = (uint8_t*)calloc(buf_size, 1);
	while ((l = knet_read(fp_remote, buf, buf_size)) != 0)
		fwrite(buf, 1, l, fp);
	free(buf);
	fclose(fp);
	knet_close(fp_remote);
}
#else
static void download_from_remote(const char *url)
{
	return;
}
#endif

bam_index_t *bam_index_load(const char *fn)
{
	bam_index_t *idx;
	idx = bam_index_load_local(fn);
	if (idx == 0 && (strstr(fn, "ftp://") == fn || strstr(fn, "http://") == fn)) {
		char *fnidx = calloc(strlen(fn) + 5, 1);
		strcat(strcpy(fnidx, fn), ".bai");
		fprintf(stderr, "[bam_index_load] attempting to download the remote index file.\n");
		download_from_remote(fnidx);
		idx = bam_index_load_local(fn);
	}
	if (idx == 0) fprintf(stderr, "[bam_index_load] fail to load BAM index.\n");
	return idx;
}

int bam_index_build2(const char *fn, const char *_fnidx)
{
	char *fnidx;
	FILE *fpidx;
	bamFile fp;
	bam_index_t *idx;
	if ((fp = bam_open(fn, "r")) == 0) {
		fprintf(stderr, "[bam_index_build2] fail to open the BAM file.\n");
		return -1;
	}
	idx = bam_index_core(fp);
	bam_close(fp);
	if (_fnidx == 0) {
		fnidx = (char*)calloc(strlen(fn) + 5, 1);
		strcpy(fnidx, fn); strcat(fnidx, ".bai");
	} else fnidx = strdup(_fnidx);
	fpidx = fopen(fnidx, "w");
	if (fpidx == 0) {
		fprintf(stderr, "[bam_index_build2] fail to create the index file.\n");
		free(fnidx);
		return -1;
	}
	bam_index_save(idx, fpidx);
	bam_index_destroy(idx);
	fclose(fpidx);
	free(fnidx);
	return 0;
}

int bam_index_build(const char *fn)
{
	return bam_index_build2(fn, 0);
}

int bam_index(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "Usage: samtools index <in.bam> [<out.index>]\n");
		return 1;
	}
	if (argc >= 3) bam_index_build2(argv[1], argv[2]);
	else bam_index_build(argv[1]);
	return 0;
}

#define MAX_BIN 37450 // =(8^6-1)/7+1

static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[MAX_BIN])
{
	int i = 0, k;
	--end;
	list[i++] = 0;
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
	for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	return i;
}

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
	uint32_t rbeg = b->core.pos;
	uint32_t rend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
	return (rend > beg && rbeg < end);
}

// bam_fetch helper function retrieves 
pair64_t * get_chunk_coordinates(const bam_index_t *idx, int tid, int beg, int end, int* cnt_off)
{
	uint16_t *bins;
	int i, n_bins, n_off;
	pair64_t *off;
	khint_t k;
	khash_t(i) *index;
	uint64_t min_off;

	bins = (uint16_t*)calloc(MAX_BIN, 2);
	n_bins = reg2bins(beg, end, bins);
	index = idx->index[tid];
	min_off = (beg>>BAM_LIDX_SHIFT >= idx->index2[tid].n)? 0 : idx->index2[tid].offset[beg>>BAM_LIDX_SHIFT];
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index))
			n_off += kh_value(index, k).n;
	}
	if (n_off == 0) {
		free(bins); return 0;
	}
	off = (pair64_t*)calloc(n_off, 16);
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index)) {
			int j;
			bam_binlist_t *p = &kh_value(index, k);
			for (j = 0; j < p->n; ++j)
				if (p->list[j].v > min_off) off[n_off++] = p->list[j];
		}
	}
	free(bins);
	{
		bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
		int l;
		ks_introsort(off, n_off, off);
		// resolve completely contained adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i)
			if (off[l].v < off[i].v)
				off[++l] = off[i];
		n_off = l + 1;
		// resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
		for (i = 1; i < n_off; ++i)
			if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
		{ // merge adjacent blocks
#if defined(BAM_TRUE_OFFSET) || defined(BAM_VIRTUAL_OFFSET16)
			for (i = 1, l = 0; i < n_off; ++i) {
#ifdef BAM_TRUE_OFFSET
				if (off[l].v + BAM_MIN_CHUNK_GAP > off[i].u) off[l].v = off[i].v;
#else
				if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
#endif
				else off[++l] = off[i];
			}
			n_off = l + 1;
#endif
		}
		bam_destroy1(b);
	}
	*cnt_off = n_off;
	return off;
}

int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
{
	int n_off;
	pair64_t *off = get_chunk_coordinates(idx, tid, beg, end, &n_off);
	if (off == 0) return 0;
	{
		// retrive alignments
		uint64_t curr_off;
		int i, ret, n_seeks;
		n_seeks = 0; i = -1; curr_off = 0;
		bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
		for (;;) {
			if (curr_off == 0 || curr_off >= off[i].v) { // then jump to the next chunk
				if (i == n_off - 1) break; // no more chunks
				if (i >= 0) assert(curr_off == off[i].v); // otherwise bug
				if (i < 0 || off[i].v != off[i+1].u) { // not adjacent chunks; then seek
					bam_seek(fp, off[i+1].u, SEEK_SET);
					curr_off = bam_tell(fp);
					++n_seeks;
				}
				++i;
			}
			if ((ret = bam_read1(fp, b)) > 0) {
				curr_off = bam_tell(fp);
				if (b->core.tid != tid || b->core.pos >= end) break; // no need to proceed
				else if (is_overlap(beg, end, b)) func(b, data);
			} else break; // end of file
		}
//		fprintf(stderr, "[bam_fetch] # seek calls: %d\n", n_seeks);
		bam_destroy1(b);
	}
	free(off);
	return 0;
}
