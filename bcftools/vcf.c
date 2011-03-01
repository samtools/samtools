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
	int max_ref;
} vcf_t;

bcf_hdr_t *vcf_hdr_read(bcf_t *bp)
{
	kstring_t meta, smpl;
	int dret;
	vcf_t *v;
	bcf_hdr_t *h;
	if (!bp->is_vcf) return bcf_hdr_read(bp);
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
			ks_tokaux_t aux;
			char *p;
			for (p = kstrtok(v->line.s, "\t\n", &aux), k = 0; p; p = kstrtok(0, 0, &aux), ++k) {
				if (k >= 9) {
					kputsn(p, aux.p - p, &smpl);
					kputc('\0', &smpl);
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
	if (strchr(mode, 'b')) return bcf_open(fn, mode);
	bp = calloc(1, sizeof(bcf_t));
	v = calloc(1, sizeof(vcf_t));
	bp->is_vcf = 1;
	bp->v = v;
	v->refhash = bcf_str2id_init();
	if (strchr(mode, 'r')) {
		v->fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
		v->ks = ks_init(v->fp);
	} else if (strchr(mode, 'w'))
		v->fpout = strcmp(fn, "-")? fopen(fn, "w") : stdout;
	return bp;
}

int vcf_dictread(bcf_t *bp, bcf_hdr_t *h, const char *fn)
{
	vcf_t *v;
	gzFile fp;
	kstream_t *ks;
	kstring_t s, rn;
	int dret;
	if (bp == 0) return -1;
	if (!bp->is_vcf) return 0;
	s.l = s.m = 0; s.s = 0;
	rn.m = rn.l = h->l_nm; rn.s = h->name;
	v = (vcf_t*)bp->v;
	fp = gzopen(fn, "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, &s, &dret) >= 0) {
		bcf_str2id_add(v->refhash, strdup(s.s));
		kputs(s.s, &rn); kputc('\0', &rn);
		if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
	}
	ks_destroy(ks);
	gzclose(fp);
	h->l_nm = rn.l; h->name = rn.s;
	bcf_hdr_sync(h);
	free(s.s);
	return 0;
}

int vcf_close(bcf_t *bp)
{
	vcf_t *v;
	if (bp == 0) return -1;
	if (!bp->is_vcf) return bcf_close(bp);
	v = (vcf_t*)bp->v;
	if (v->fp) {
		ks_destroy(v->ks);
		gzclose(v->fp);
	}
	if (v->fpout) fclose(v->fpout);
	free(v->line.s);
	bcf_str2id_thorough_destroy(v->refhash);
	free(v);
	free(bp);
	return 0;
}

int vcf_hdr_write(bcf_t *bp, const bcf_hdr_t *h)
{
	vcf_t *v = (vcf_t*)bp->v;
	int i, has_ver = 0;
	if (!bp->is_vcf) return bcf_hdr_write(bp, h);
	if (h->l_txt > 0) {
		if (strstr(h->txt, "##fileformat=")) has_ver = 1;
		if (has_ver == 0) fprintf(v->fpout, "##fileformat=VCFv4.1\n");
		fwrite(h->txt, 1, h->l_txt - 1, v->fpout);
	}
	if (h->l_txt == 0) fprintf(v->fpout, "##fileformat=VCFv4.1\n");
	fprintf(v->fpout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (i = 0; i < h->n_smpl; ++i)
		fprintf(v->fpout, "\t%s", h->sns[i]);
	fputc('\n', v->fpout);
	return 0;
}

int vcf_write(bcf_t *bp, bcf_hdr_t *h, bcf1_t *b)
{
	vcf_t *v = (vcf_t*)bp->v;
	extern void bcf_fmt_core(const bcf_hdr_t *h, bcf1_t *b, kstring_t *s);
	if (!bp->is_vcf) return bcf_write(bp, h, b);
	bcf_fmt_core(h, b, &v->line);
	fwrite(v->line.s, 1, v->line.l, v->fpout);
	fputc('\n', v->fpout);
	return v->line.l + 1;
}

int vcf_read(bcf_t *bp, bcf_hdr_t *h, bcf1_t *b)
{
	int dret, k, i, sync = 0;
	vcf_t *v = (vcf_t*)bp->v;
	char *p, *q;
	kstring_t str, rn;
	ks_tokaux_t aux, a2;
	if (!bp->is_vcf) return bcf_read(bp, h, b);
	v->line.l = 0;
	str.l = 0; str.m = b->m_str; str.s = b->str;
	rn.l = rn.m = h->l_nm; rn.s = h->name;
	if (ks_getuntil(v->ks, '\n', &v->line, &dret) < 0) return -1;
	b->n_smpl = h->n_smpl;
	for (p = kstrtok(v->line.s, "\t", &aux), k = 0; p; p = kstrtok(0, 0, &aux), ++k) {
		*(char*)aux.p = 0;
		if (k == 0) { // ref
			int tid = bcf_str2id(v->refhash, p);
			if (tid < 0) {
				tid = bcf_str2id_add(v->refhash, strdup(p));
				kputs(p, &rn); kputc('\0', &rn);
				sync = 1;
			}
			b->tid = tid;
		} else if (k == 1) { // pos
			b->pos = atoi(p) - 1;
		} else if (k == 5) { // qual
			b->qual = (p[0] >= '0' && p[0] <= '9')? atof(p) : 0;
		} else if (k <= 8) { // variable length strings
			kputs(p, &str); kputc('\0', &str);
			b->l_str = str.l; b->m_str = str.m; b->str = str.s;
			if (k == 8) bcf_sync(b);
		} else { // k > 9
			if (strncmp(p, "./.", 3) == 0) {
				for (i = 0; i < b->n_gi; ++i) {
					if (b->gi[i].fmt == bcf_str2int("GT", 2)) {
						((uint8_t*)b->gi[i].data)[k-9] = 1<<7;
					} else if (b->gi[i].fmt == bcf_str2int("GQ", 2)) {
						((uint8_t*)b->gi[i].data)[k-9] = 0;
					} else if (b->gi[i].fmt == bcf_str2int("SP", 2)) {
						((int32_t*)b->gi[i].data)[k-9] = 0;
					} else if (b->gi[i].fmt == bcf_str2int("DP", 2)) {
						((uint16_t*)b->gi[i].data)[k-9] = 0;
					} else if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
						int y = b->n_alleles * (b->n_alleles + 1) / 2;
						memset((uint8_t*)b->gi[i].data + (k - 9) * y, 0, y);
					} else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
						int y = b->n_alleles * (b->n_alleles + 1) / 2;
						memset((float*)b->gi[i].data + (k - 9) * y, 0, y * 4);
					}
				}
				goto endblock;
			}
			for (q = kstrtok(p, ":", &a2), i = 0; q && i < b->n_gi; q = kstrtok(0, 0, &a2), ++i) {
				if (b->gi[i].fmt == bcf_str2int("GT", 2)) {
					((uint8_t*)b->gi[i].data)[k-9] = (q[0] - '0')<<3 | (q[2] - '0') | (q[1] == '/'? 0 : 1) << 6;
				} else if (b->gi[i].fmt == bcf_str2int("GQ", 2)) {
					double _x = strtod(q, &q);
					int x = (int)(_x + .499);
					if (x > 255) x = 255;
					((uint8_t*)b->gi[i].data)[k-9] = x;
				} else if (b->gi[i].fmt == bcf_str2int("SP", 2)) {
					int x = strtol(q, &q, 10);
					if (x > 0xffff) x = 0xffff;
					((uint32_t*)b->gi[i].data)[k-9] = x;
				} else if (b->gi[i].fmt == bcf_str2int("DP", 2)) {
					int x = strtol(q, &q, 10);
					if (x > 0xffff) x = 0xffff;
					((uint16_t*)b->gi[i].data)[k-9] = x;
				} else if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
					int x, y, j;
					uint8_t *data = (uint8_t*)b->gi[i].data;
					y = b->n_alleles * (b->n_alleles + 1) / 2;
					for (j = 0; j < y; ++j) {
						x = strtol(q, &q, 10);
						if (x > 255) x = 255;
						data[(k-9) * y + j] = x;
						++q;
					}
				} else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
					int j, y;
					float x, *data = (float*)b->gi[i].data;
					y = b->n_alleles * (b->n_alleles + 1) / 2;
					for (j = 0; j < y; ++j) {
						x = strtod(q, &q);
						data[(k-9) * y + j] = x > 0? -x/10. : x;
						++q;
					}
				}
			}
		endblock: i = i;
		}
	}
	h->l_nm = rn.l; h->name = rn.s;
	if (sync) bcf_hdr_sync(h);
	return v->line.l + 1;
}
