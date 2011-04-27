#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include "bcf.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

int bcfview(int argc, char *argv[]);
int bcf_main_index(int argc, char *argv[]);

#define BUF_SIZE 0x10000

int bcf_cat(int n, char * const *fn)
{
	int i;
	bcf_t *out;
	uint8_t *buf;
	buf = malloc(BUF_SIZE);
	out = bcf_open("-", "w");
	for (i = 0; i < n; ++i) {
		bcf_t *in;
		bcf_hdr_t *h;
		off_t end;
		struct stat s;
		in = bcf_open(fn[i], "r");
		h = bcf_hdr_read(in);
		if (i == 0) bcf_hdr_write(out, h);
		bcf_hdr_destroy(h);
#ifdef _USE_KNETFILE
		fstat(knet_fileno(in->fp->x.fpr), &s);
		end = s.st_size - 28;
		while (knet_tell(in->fp->x.fpr) < end) {
			int size = knet_tell(in->fp->x.fpr) + BUF_SIZE < end? BUF_SIZE : end - knet_tell(in->fp->x.fpr);
			knet_read(in->fp->x.fpr, buf, size);
			fwrite(buf, 1, size, out->fp->x.fpw);
		}
#else
		abort(); // FIXME: not implemented
#endif
		bcf_close(in);
	}
	bcf_close(out);
	free(buf);
	return 0;
}

extern double bcf_pair_freq(const bcf1_t *b0, const bcf1_t *b1, double f[4]);

int bcf_main_ldpair(int argc, char *argv[])
{
	bcf_t *fp;
	bcf_hdr_t *h;
	bcf1_t *b0, *b1;
	bcf_idx_t *idx;
	kstring_t str;
	void *str2id;
	gzFile fplist;
	kstream_t *ks;
	int dret, lineno = 0;
	if (argc < 3) {
		fprintf(stderr, "Usage: bcftools ldpair <in.bcf> <in.list>\n");
		return 1;
	}
	fplist = gzopen(argv[2], "rb");
	ks = ks_init(fplist);
	memset(&str, 0, sizeof(kstring_t));
	fp = bcf_open(argv[1], "rb");
	h = bcf_hdr_read(fp);
	str2id = bcf_build_refhash(h);
	idx = bcf_idx_load(argv[1]);
	if (idx == 0) {
		fprintf(stderr, "[%s] No bcf index is found. Abort!\n", __func__);
		return 1;
	}
	b0 = calloc(1, sizeof(bcf1_t));
	b1 = calloc(1, sizeof(bcf1_t));
	while (ks_getuntil(ks, '\n', &str, &dret) >= 0) {
		char *p, *q;
		int k;
		int tid0 = -1, tid1 = -1, pos0 = -1, pos1 = -1;
		++lineno;
		for (p = q = str.s, k = 0; *p; ++p) {
			if (*p == ' ' || *p == '\t') {
				*p = '\0';
				if (k == 0) tid0 = bcf_str2id(str2id, q);
				else if (k == 1) pos0 = atoi(q) - 1;
				else if (k == 2) tid1 = strcmp(q, "=")? bcf_str2id(str2id, q) : tid0;
				else if (k == 3) pos1 = atoi(q) - 1;
				q = p + 1;
				++k;
			}
		}
		if (k == 3) pos1 = atoi(q) - 1;
		if (tid0 >= 0 && tid1 >= 0 && pos0 >= 0 && pos1 >= 0) {
			uint64_t off;
			double r, f[4];
			off = bcf_idx_query(idx, tid0, pos0);
			bgzf_seek(fp->fp, off, SEEK_SET);
			while (bcf_read(fp, h, b0) >= 0 && b0->pos != pos0);
			off = bcf_idx_query(idx, tid1, pos1);
			bgzf_seek(fp->fp, off, SEEK_SET);
			while (bcf_read(fp, h, b1) >= 0 && b1->pos != pos1);
			r = bcf_pair_freq(b0, b1, f);
			r *= r;
			printf("%s\t%d\t%s\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n", h->ns[tid0], pos0+1, h->ns[tid1], pos1+1,
				r, f[0], f[1], f[2], f[3]);
		} //else fprintf(stderr, "[%s] Parse error at line %d.\n", __func__, lineno);
	}
	bcf_destroy(b0); bcf_destroy(b1);
	bcf_idx_destroy(idx);
	bcf_str2id_destroy(str2id);
	bcf_hdr_destroy(h);
	bcf_close(fp);
	free(str.s);
	ks_destroy(ks);
	gzclose(fplist);
	return 0;
}

int bcf_main_ld(int argc, char *argv[])
{
	bcf_t *fp;
	bcf_hdr_t *h;
	bcf1_t **b, *b0;
	int i, j, m, n;
	double f[4];
	if (argc == 1) {
		fprintf(stderr, "Usage: bcftools ld <in.bcf>\n");
		return 1;
	}
	fp = bcf_open(argv[1], "rb");
	h = bcf_hdr_read(fp);
	// read the entire BCF
	m = n = 0; b = 0;
	b0 = calloc(1, sizeof(bcf1_t));
	while (bcf_read(fp, h, b0) >= 0) {
		if (m == n) {
			m = m? m<<1 : 16;
			b = realloc(b, sizeof(void*) * m);
		}
		b[n] = calloc(1, sizeof(bcf1_t));
		bcf_cpy(b[n++], b0);
	}
	bcf_destroy(b0);
	// compute pair-wise r^2
	printf("%d\n", n); // the number of loci
	for (i = 0; i < n; ++i) {
		printf("%s:%d", h->ns[b[i]->tid], b[i]->pos + 1);
		for (j = 0; j < i; ++j) {
			double r = bcf_pair_freq(b[i], b[j], f);
			printf("\t%.3f", r*r);
		}
		printf("\t1.000\n");
	}
	// free
	for (i = 0; i < n; ++i) bcf_destroy(b[i]);
	free(b);
	bcf_hdr_destroy(h);
	bcf_close(fp);
	return 0;
}

int main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: bcftools (Tools for data in the VCF/BCF formats)\n");
		fprintf(stderr, "Version: %s\n\n", BCF_VERSION);
		fprintf(stderr, "Usage:   bcftools <command> <arguments>\n\n");
		fprintf(stderr, "Command: view      print, extract, convert and call SNPs from BCF\n");
		fprintf(stderr, "         index     index BCF\n");
		fprintf(stderr, "         cat       concatenate BCFs\n");
		fprintf(stderr, "         ld        compute all-pair r^2\n");
		fprintf(stderr, "         ldpair    compute r^2 between requested pairs\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "view") == 0) return bcfview(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) return bcf_main_index(argc-1, argv+1);
	else if (strcmp(argv[1], "ld") == 0) return bcf_main_ld(argc-1, argv+1);
	else if (strcmp(argv[1], "ldpair") == 0) return bcf_main_ldpair(argc-1, argv+1);
	else if (strcmp(argv[1], "cat") == 0) return bcf_cat(argc-2, argv+2); // cat is different ...
	else {
		fprintf(stderr, "[main] Unrecognized command.\n");
		return 1;
	}
	return 0;
}
