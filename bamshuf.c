#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sam.h"
#include "ksort.h"

#define DEF_CLEVEL 1

static inline unsigned hash_Wang(unsigned key)
{
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}

static inline unsigned hash_X31_Wang(const char *s)
{
	unsigned h = *s;
	if (h) {
		for (++s ; *s; ++s) h = (h << 5) - h + *s;
		return hash_Wang(h);
	} else return 0;
}

typedef struct {
	unsigned key;
	bam1_t *b;
} elem_t;

static inline int elem_lt(elem_t x, elem_t y)
{
	if (x.key < y.key) return 1;
	if (x.key == y.key) {
		int t;
		t = strcmp(bam_get_qname(x.b), bam_get_qname(y.b));
		if (t < 0) return 1;
		return (t == 0 && ((x.b->core.flag>>6&3) < (y.b->core.flag>>6&3)));
	} else return 0;
}

KSORT_INIT(bamshuf, elem_t, elem_lt)

static void bamshuf(const char *fn, int n_files, const char *pre, int clevel, int is_stdout)
{
	BGZF *fp, *fpw, **fpt;
	char **fnt, modew[8];
	bam1_t *b;
	int i, l;
	bam_hdr_t *h;
	int64_t *cnt;

	// split
	fp = strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
	assert(fp);
	h = bam_hdr_read(fp);
	fnt = (char**)calloc(n_files, sizeof(void*));
	fpt = (BGZF**)calloc(n_files, sizeof(void*));
	cnt = (int64_t*)calloc(n_files, 8);
	l = strlen(pre);
	for (i = 0; i < n_files; ++i) {
		fnt[i] = (char*)calloc(l + 10, 1);
		sprintf(fnt[i], "%s.%.4d.bam", pre, i);
		fpt[i] = bgzf_open(fnt[i], "w1");
		bam_hdr_write(fpt[i], h);
	}
	b = bam_init1();
	while (bam_read1(fp, b) >= 0) {
		uint32_t x;
		x = hash_X31_Wang(bam_get_qname(b)) % n_files;
		bam_write1(fpt[x], b);
		++cnt[x];
	}
	bam_destroy1(b);
	for (i = 0; i < n_files; ++i) bgzf_close(fpt[i]);
	free(fpt);
	bgzf_close(fp);
	// merge
	sprintf(modew, "w%d", (clevel >= 0 && clevel <= 9)? clevel : DEF_CLEVEL);
	if (!is_stdout) { // output to a file
		char *fnw = (char*)calloc(l + 5, 1);
		sprintf(fnw, "%s.bam", pre);
		fpw = bgzf_open(fnw, modew);
		free(fnw);
	} else fpw = bgzf_dopen(fileno(stdout), modew); // output to stdout
	bam_hdr_write(fpw, h);
	bam_hdr_destroy(h);
	for (i = 0; i < n_files; ++i) {
		int64_t j, c = cnt[i];
		elem_t *a;
		fp = bgzf_open(fnt[i], "r");
		bam_hdr_destroy(bam_hdr_read(fp));
		a = (elem_t*)calloc(c, sizeof(elem_t));
		for (j = 0; j < c; ++j) {
			a[j].b = bam_init1();
			assert(bam_read1(fp, a[j].b) >= 0);
			a[j].key = hash_X31_Wang(bam_get_qname(a[j].b));
		}
		bgzf_close(fp);
		unlink(fnt[i]);
		free(fnt[i]);
		ks_introsort(bamshuf, c, a);
		for (j = 0; j < c; ++j) {
			bam_write1(fpw, a[j].b);
			bam_destroy1(a[j].b);
		}
		free(a);
	}
	bgzf_close(fpw);
	free(fnt); free(cnt);
}

int main_bamshuf(int argc, char *argv[])
{
	int c, n_files = 64, clevel = DEF_CLEVEL, is_stdout = 0, is_un = 0;
	while ((c = getopt(argc, argv, "n:l:uO")) >= 0) {
		switch (c) {
		case 'n': n_files = atoi(optarg); break;
		case 'l': clevel = atoi(optarg); break;
		case 'u': is_un = 1; break;
		case 'O': is_stdout = 1; break;
		}
	}
	if (is_un) clevel = 0;
	if (optind + 2 > argc) {
		fprintf(stderr, "\nUsage:   bamshuf [-Ou] [-n nFiles] [-c cLevel] <in.bam> <out.prefix>\n\n");
		fprintf(stderr, "Options: -O      output to stdout\n");
		fprintf(stderr, "         -u      uncompressed BAM output\n");
		fprintf(stderr, "         -l INT  compression level [%d]\n", DEF_CLEVEL);
		fprintf(stderr, "         -n INT  number of temporary files [%d]\n", n_files);
		fprintf(stderr, "\n");
		return 1;
	}
	bamshuf(argv[optind], n_files, argv[optind+1], clevel, is_stdout);
	return 0;
}
