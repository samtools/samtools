#include <assert.h>
#include <ctype.h>
#include <sys/stat.h>
#include "bam_endian.h"
#include "kstring.h"
#include "bcf.h"
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

#define TAD_LIDX_SHIFT 13

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} bcf_lidx_t;

struct __bcf_idx_t {
	int32_t n;
	bcf_lidx_t *index2;
};

/************
 * indexing *
 ************/

static inline void insert_offset2(bcf_lidx_t *index2, int _beg, int _end, uint64_t offset)
{
	int i, beg, end;
	beg = _beg >> TAD_LIDX_SHIFT;
	end = (_end - 1) >> TAD_LIDX_SHIFT;
	if (index2->m < end + 1) {
		int old_m = index2->m;
		index2->m = end + 1;
		kroundup32(index2->m);
		index2->offset = (uint64_t*)realloc(index2->offset, index2->m * 8);
		memset(index2->offset + old_m, 0, 8 * (index2->m - old_m));
	}
	if (beg == end) {
		if (index2->offset[beg] == 0) index2->offset[beg] = offset;
	} else {
		for (i = beg; i <= end; ++i)
			if (index2->offset[i] == 0) index2->offset[i] = offset;
	}
	if (index2->n < end + 1) index2->n = end + 1;
}

bcf_idx_t *bcf_idx_core(bcf_t *bp, bcf_hdr_t *h)
{
	bcf_idx_t *idx;
	int32_t last_coor, last_tid;
	uint64_t last_off;
	kstring_t *str;
	BGZF *fp = bp->fp;
	bcf1_t *b;
	int ret;

	b = calloc(1, sizeof(bcf1_t));
	str = calloc(1, sizeof(kstring_t));
	idx = (bcf_idx_t*)calloc(1, sizeof(bcf_idx_t));
	idx->n = h->n_ref;
	idx->index2 = calloc(h->n_ref, sizeof(bcf_lidx_t));

	last_tid = 0xffffffffu;
	last_off = bgzf_tell(fp); last_coor = 0xffffffffu;
	while ((ret = bcf_read(bp, h, b)) > 0) {
		int end, tmp;
		if (last_tid != b->tid) { // change of chromosomes
			last_tid = b->tid;
		} else if (last_coor > b->pos) {
			fprintf(stderr, "[bcf_idx_core] the input is out of order\n");
			free(str->s); free(str); free(idx); bcf_destroy(b);
			return 0;
		}
		tmp = strlen(b->ref);
		end = b->pos + (tmp > 0? tmp : 1);
		insert_offset2(&idx->index2[b->tid], b->pos, end, last_off);
		last_off = bgzf_tell(fp);
		last_coor = b->pos;
	}
	free(str->s); free(str); bcf_destroy(b);
	return idx;
}

void bcf_idx_destroy(bcf_idx_t *idx)
{
	int i;
	if (idx == 0) return;
	for (i = 0; i < idx->n; ++i) free(idx->index2[i].offset);
	free(idx->index2);
	free(idx);
}

/******************
 * index file I/O *
 ******************/

void bcf_idx_save(const bcf_idx_t *idx, BGZF *fp)
{
	int32_t i, ti_is_be;
	ti_is_be = bam_is_big_endian();
	bgzf_write(fp, "BCI\4", 4);
	if (ti_is_be) {
		uint32_t x = idx->n;
		bgzf_write(fp, bam_swap_endian_4p(&x), 4);
	} else bgzf_write(fp, &idx->n, 4);
	for (i = 0; i < idx->n; ++i) {
		bcf_lidx_t *index2 = idx->index2 + i;
		// write linear index (index2)
		if (ti_is_be) {
			int x = index2->n;
			bgzf_write(fp, bam_swap_endian_4p(&x), 4);
		} else bgzf_write(fp, &index2->n, 4);
		if (ti_is_be) { // big endian
			int x;
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
			bgzf_write(fp, index2->offset, 8 * index2->n);
			for (x = 0; (int)x < index2->n; ++x)
				bam_swap_endian_8p(&index2->offset[x]);
		} else bgzf_write(fp, index2->offset, 8 * index2->n);
	}
}

static bcf_idx_t *bcf_idx_load_core(BGZF *fp)
{
	int i, ti_is_be;
	char magic[4];
	bcf_idx_t *idx;
	ti_is_be = bam_is_big_endian();
	if (fp == 0) {
		fprintf(stderr, "[%s] fail to load index.\n", __func__);
		return 0;
	}
	bgzf_read(fp, magic, 4);
	if (strncmp(magic, "BCI\4", 4)) {
		fprintf(stderr, "[%s] wrong magic number.\n", __func__);
		return 0;
	}
	idx = (bcf_idx_t*)calloc(1, sizeof(bcf_idx_t));	
	bgzf_read(fp, &idx->n, 4);
	if (ti_is_be) bam_swap_endian_4p(&idx->n);
	idx->index2 = (bcf_lidx_t*)calloc(idx->n, sizeof(bcf_lidx_t));
	for (i = 0; i < idx->n; ++i) {
		bcf_lidx_t *index2 = idx->index2 + i;
		int j;
		bgzf_read(fp, &index2->n, 4);
		if (ti_is_be) bam_swap_endian_4p(&index2->n);
		index2->m = index2->n;
		index2->offset = (uint64_t*)calloc(index2->m, 8);
		bgzf_read(fp, index2->offset, index2->n * 8);
		if (ti_is_be)
			for (j = 0; j < index2->n; ++j) bam_swap_endian_8p(&index2->offset[j]);
	}
	return idx;
}

bcf_idx_t *bcf_idx_load_local(const char *fnidx)
{
	BGZF *fp;
	fp = bgzf_open(fnidx, "r");
	if (fp) {
		bcf_idx_t *idx = bcf_idx_load_core(fp);
		bgzf_close(fp);
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

static char *get_local_version(const char *fn)
{
    struct stat sbuf;
	char *fnidx = (char*)calloc(strlen(fn) + 5, 1);
	strcat(strcpy(fnidx, fn), ".bci");
	if ((strstr(fnidx, "ftp://") == fnidx || strstr(fnidx, "http://") == fnidx)) {
		char *p, *url;
		int l = strlen(fnidx);
		for (p = fnidx + l - 1; p >= fnidx; --p)
			if (*p == '/') break;
		url = fnidx; fnidx = strdup(p + 1);
		if (stat(fnidx, &sbuf) == 0) {
			free(url);
			return fnidx;
		}
		fprintf(stderr, "[%s] downloading the index file...\n", __func__);
		download_from_remote(url);
		free(url);
	}
    if (stat(fnidx, &sbuf) == 0) return fnidx;
	free(fnidx); return 0;
}

bcf_idx_t *bcf_idx_load(const char *fn)
{
	bcf_idx_t *idx;
    char *fname = get_local_version(fn);
	if (fname == 0) return 0;
	idx = bcf_idx_load_local(fname);
    free(fname);
	return idx;
}

int bcf_idx_build2(const char *fn, const char *_fnidx)
{
	char *fnidx;
	BGZF *fpidx;
	bcf_t *bp;
	bcf_idx_t *idx;
	bcf_hdr_t *h;
	if ((bp = bcf_open(fn, "r")) == 0) {
		fprintf(stderr, "[bcf_idx_build2] fail to open the BAM file.\n");
		return -1;
	}
	h = bcf_hdr_read(bp);
	idx = bcf_idx_core(bp, h);
	bcf_close(bp);
	if (_fnidx == 0) {
		fnidx = (char*)calloc(strlen(fn) + 5, 1);
		strcpy(fnidx, fn); strcat(fnidx, ".bci");
	} else fnidx = strdup(_fnidx);
	fpidx = bgzf_open(fnidx, "w");
	if (fpidx == 0) {
		fprintf(stderr, "[bcf_idx_build2] fail to create the index file.\n");
		free(fnidx);
		return -1;
	}
	bcf_idx_save(idx, fpidx);
	bcf_idx_destroy(idx);
	bgzf_close(fpidx);
	free(fnidx);
	return 0;
}

int bcf_idx_build(const char *fn)
{
	return bcf_idx_build2(fn, 0);
}

/********************************************
 * parse a region in the format chr:beg-end *
 ********************************************/

int bcf_parse_region(void *str2id, const char *str, int *tid, int *begin, int *end)
{
	char *s, *p;
	int i, l, k;
	l = strlen(str);
	p = s = (char*)malloc(l+1);
	/* squeeze out "," */
	for (i = k = 0; i != l; ++i)
		if (str[i] != ',' && !isspace(str[i])) s[k++] = str[i];
	s[k] = 0;
	for (i = 0; i != k; ++i) if (s[i] == ':') break;
	s[i] = 0;
	if ((*tid = bcf_str2id(str2id, s)) < 0) {
		free(s);
		return -1;
	}
	if (i == k) { /* dump the whole sequence */
		*begin = 0; *end = 1<<29; free(s);
		return 0;
	}
	for (p = s + i + 1; i != k; ++i) if (s[i] == '-') break;
	*begin = atoi(p);
	if (i < k) {
		p = s + i + 1;
		*end = atoi(p);
	} else *end = 1<<29;
	if (*begin > 0) --*begin;
	free(s);
	if (*begin > *end) return -1;
	return 0;
}

/*******************************
 * retrieve a specified region *
 *******************************/

uint64_t bcf_idx_query(const bcf_idx_t *idx, int tid, int beg)
{
	uint64_t min_off, *offset;
	int i;
	if (beg < 0) beg = 0;
	offset = idx->index2[tid].offset;
	for (i = beg>>TAD_LIDX_SHIFT; i < idx->index2[tid].n && offset[i] == 0; ++i);
	min_off = (i == idx->index2[tid].n)? offset[idx->index2[tid].n-1] : offset[i];
	return min_off;
}

int bcf_main_index(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "Usage: bcftools index <in.bcf>\n");
		return 1;
	}
	bcf_idx_build(argv[1]);
	return 0;
}
