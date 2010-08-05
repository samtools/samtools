#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bcf.h"
#include "kstring.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 4096)

typedef struct {
	gzFile fp;
	FILE *fpout;
	kstream_t *ks;
	void *refhash;
	kstring_t line;
} vcf_t;

bcf_hdr_t *vcf_hdr_read(bcf_t *bp)
{
	kstring_t meta, smpl;
	int dret;
	vcf_t *v;
	bcf_hdr_t *h;
	if (!bp->is_vcf) return 0;
	h = calloc(1, sizeof(bcf_hdr_t));
	v = (vcf_t*)bp->v;
	v->line.l = 0;
	memset(&meta, 0, sizeof(kstring_t));
	memset(&smpl, 0, sizeof(kstring_t));
	while (ks_getuntil(v->ks, '\n', &v->line, &dret) >= 0) {
		if (v->line.l < 2) continue;
		if (v->line.s[0] != '#') return 0; // no sample line
		if (v->line.s[0] == '#' && v->line.s[1] == '#') {
			kputsn(v->line.s, v->line.l, &meta); kputc('\n', &meta);
		} else if (v->line.s[0] == '#') {
			int k;
			char *p, *q, *r;
			for (q = v->line.s, p = q + 1, k = 0; *p; ++p) {
				if (*p == '\t' || *(p+1) == 0) {
					r = *(p+1) == 0? p+1 : p;
					if (k >= 9) {
						kputsn(q, r - q, &smpl);
						kputc('\0', &smpl);
					}
					q = p + 1; ++k;
				}
			}
			break;
		}
	}
	kputc('\0', &meta);
	h->name = 0;
	h->sname = smpl.s; h->l_smpl = smpl.l;
	h->txt = meta.s; h->l_txt = meta.l;
	bcf_hdr_sync(h);
	return h;
}

bcf_t *vcf_open(const char *fn, const char *mode)
{
	bcf_t *bp;
	vcf_t *v;
	bp = calloc(1, sizeof(bcf_t));
	v = calloc(1, sizeof(vcf_t));
	bp->is_vcf = 1;
	bp->v = v;
	if (strchr(mode, 'r')) {
		v->fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
		v->ks = ks_init(v->fp);
	} else if (strchr(mode, 'w'))
		v->fpout = strcmp(fn, "-")? fopen(fn, "w") : fdopen(fileno(stdout), "w");
	return bp;
}

void bcf_hdr_clear(bcf_hdr_t *b);

int vcf_close(bcf_t *bp)
{
	vcf_t *v;
	if (bp == 0) return -1;
	if (bp->v == 0) return -1;
	v = (vcf_t*)bp->v;
	if (v->fp) {
		ks_destroy(v->ks);
		gzclose(v->fp);
	}
	if (v->fpout) fclose(v->fpout);
	free(v->line.s);
	free(v);
	free(bp);
	return 0;
}

int vcf_hdr_write(bcf_t *bp, const bcf_hdr_t *h)
{
	vcf_t *v = (vcf_t*)bp->v;
	int i;
	if (v == 0 || v->fpout == 0) return -1;
	fwrite(h->txt, 1, h->l_txt, v->fpout);
	fprintf(v->fpout, "#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (i = 0; i < h->n_smpl; ++i)
		fprintf(v->fpout, "\t%s", h->sns[i]);
	fputc('\n', v->fpout);
	return 0;
}

int vcf_read(bcf_t *bp, bcf1_t *b)
{
	int dret;
	vcf_t *v = (vcf_t*)bp->v;
	v->line.l = 0;
	if (ks_getuntil(v->ks, '\n', &v->line, &dret) < 0) return -1;
	return v->line.l + 1;
}

int vcf_test(int argc, char *argv[])
{
	bcf_t *bp, *bpout;
	bcf_hdr_t *h;
	bp = vcf_open(argv[1], "r");
	bpout = vcf_open("-", "w");
	h = vcf_hdr_read(bpout);
	vcf_hdr_write(bpout, h);
	vcf_close(bp);
	return 0;
}
