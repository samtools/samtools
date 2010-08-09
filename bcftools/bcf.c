#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "kstring.h"
#include "bcf.h"

void bcf_hdr_clear(bcf_hdr_t *b)
{
	free(b->name); free(b->sname); free(b->txt); free(b->ns); free(b->sns);
	memset(b, 0, sizeof(bcf_hdr_t));
}

bcf_t *bcf_open(const char *fn, const char *mode)
{
	bcf_t *b;
	b = calloc(1, sizeof(bcf_t));
	if (strchr(mode, 'w')) {
		b->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_fdopen(fileno(stdout), mode);
	} else {
		b->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_fdopen(fileno(stdin), mode);
	}
	b->fp->owned_file = 1;
	return b;
}

int bcf_close(bcf_t *b)
{
	int ret;
	if (b == 0) return 0;
	ret = bgzf_close(b->fp);
	free(b);
	return ret;
}

int bcf_hdr_write(bcf_t *b, const bcf_hdr_t *h)
{
	if (b == 0 || h == 0) return -1;
	bgzf_write(b->fp, "BCF\4", 4);
	bgzf_write(b->fp, &h->l_nm, 4);
	bgzf_write(b->fp, h->name, h->l_nm);
	bgzf_write(b->fp, &h->l_smpl, 4);
	bgzf_write(b->fp, h->sname, h->l_smpl);
	bgzf_write(b->fp, &h->l_txt, 4);
	bgzf_write(b->fp, h->txt, h->l_txt);
	bgzf_flush(b->fp);
	return 16 + h->l_nm + h->l_smpl + h->l_txt;
}

bcf_hdr_t *bcf_hdr_read(bcf_t *b)
{
	uint8_t magic[4];
	bcf_hdr_t *h;
	if (b == 0) return 0;
	h = calloc(1, sizeof(bcf_hdr_t));
	bgzf_read(b->fp, magic, 4);
	bgzf_read(b->fp, &h->l_nm, 4);
	h->name = malloc(h->l_nm);
	bgzf_read(b->fp, h->name, h->l_nm);
	bgzf_read(b->fp, &h->l_smpl, 4);
	h->sname = malloc(h->l_smpl);
	bgzf_read(b->fp, h->sname, h->l_smpl);
	bgzf_read(b->fp, &h->l_txt, 4);
	h->txt = malloc(h->l_txt);
	bgzf_read(b->fp, h->txt, h->l_txt);
	bcf_hdr_sync(h);
	return h;
}

void bcf_hdr_destroy(bcf_hdr_t *h)
{
	if (h == 0) return;
	free(h->name); free(h->sname); free(h->txt); free(h->ns); free(h->sns);
	free(h);
}

static inline char **cnt_null(int l, char *str, int *_n)
{
	int n = 0;
	char *p, **list;
	*_n = 0;
	if (l == 0 || str == 0) return 0;
	for (p = str; p != str + l; ++p)
		if (*p == 0) ++n;
	*_n = n;
	list = calloc(n, sizeof(void*));
	list[0] = str;
	for (p = str, n = 1; p < str + l - 1; ++p)
		if (*p == 0) list[n++] = p + 1;
	return list;
}

int bcf_hdr_sync(bcf_hdr_t *b)
{
	if (b == 0) return -1;
	if (b->ns) free(b->ns);
	if (b->sns) free(b->sns);
	if (b->l_nm) b->ns = cnt_null(b->l_nm, b->name, &b->n_ref);
	else b->ns = 0, b->n_ref = 0;
	b->sns = cnt_null(b->l_smpl, b->sname, &b->n_smpl);
	return 0;
}

int bcf_sync(int n_smpl, bcf1_t *b)
{
	char *p, *tmp[5];
	int i, n;
	ks_tokaux_t aux;
	// set ref, alt, flt, info, fmt
	b->ref = b->alt = b->flt = b->info = b->fmt = 0;
	for (p = b->str, n = 0; p < b->str + b->l_str; ++p)
		if (*p == 0 && p+1 != b->str + b->l_str) tmp[n++] = p + 1;
	if (n != 5) return -1;
	b->ref = tmp[0]; b->alt = tmp[1]; b->flt = tmp[2]; b->info = tmp[3]; b->fmt = tmp[4];
	// set n_alleles
	for (p = b->alt, n = 1; *p; ++p)
		if (*p == ',') ++n;
	b->n_alleles = n + 1;
	// set n_gi and gi[i].fmt
	for (p = b->fmt, n = 1; *p; ++p)
		if (*p == ':') ++n;
	if (n > b->m_gi) {
		int old_m = b->m_gi;
		b->m_gi = n;
		kroundup32(b->m_gi);
		b->gi = realloc(b->gi, b->m_gi * sizeof(bcf_ginfo_t));
		memset(b->gi + old_m, 0, (b->m_gi - old_m) * sizeof(bcf_ginfo_t));
	}
	b->n_gi = n;
	for (p = kstrtok(b->fmt, ":", &aux), n = 0; p; p = kstrtok(0, 0, &aux))
		b->gi[n++].fmt = bcf_str2int(p, aux.p - p);
	// set gi[i].len
	for (i = 0; i < b->n_gi; ++i) {
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
			b->gi[i].len = b->n_alleles * (b->n_alleles + 1) / 2;
		} else if (b->gi[i].fmt == bcf_str2int("DP", 2) || b->gi[i].fmt == bcf_str2int("HQ", 2)) {
			b->gi[i].len = 2;
		} else if (b->gi[i].fmt == bcf_str2int("GQ", 2) || b->gi[i].fmt == bcf_str2int("GT", 2)) {
			b->gi[i].len = 1;
		} else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
			b->gi[i].len = 4;
		}
		b->gi[i].data = realloc(b->gi[i].data, n_smpl * b->gi[i].len);
	}
	return 0;
}

int bcf_write(bcf_t *bp, const bcf_hdr_t *h, const bcf1_t *b)
{
	uint32_t x;
	int i, l = 0;
	if (b == 0) return -1;
	bgzf_write(bp->fp, &b->tid, 4);
	bgzf_write(bp->fp, &b->pos, 4);
	x = b->qual<<24 | b->l_str;
	bgzf_write(bp->fp, &x, 4);
	bgzf_write(bp->fp, b->str, b->l_str);
	l = 12 + b->l_str;
	for (i = 0; i < b->n_gi; ++i) {
		bgzf_write(bp->fp, b->gi[i].data, b->gi[i].len * h->n_smpl);
		l += b->gi[i].len * h->n_smpl;
	}
	return l;
}

int bcf_read(bcf_t *bp, const bcf_hdr_t *h, bcf1_t *b)
{
	int i, l = 0;
	uint32_t x;
	if (b == 0) return -1;
	if (bgzf_read(bp->fp, &b->tid, 4) == 0) return -1;
	bgzf_read(bp->fp, &b->pos, 4);
	bgzf_read(bp->fp, &x, 4);
	b->qual = x >> 24; b->l_str = x << 8 >> 8;
	if (b->l_str > b->m_str) {
		b->m_str = b->l_str;
		kroundup32(b->m_str);
		b->str = realloc(b->str, b->m_str);
	}
	bgzf_read(bp->fp, b->str, b->l_str);
	l = 12 + b->l_str;
	bcf_sync(h->n_smpl, b);
	for (i = 0; i < b->n_gi; ++i) {
		bgzf_read(bp->fp, b->gi[i].data, b->gi[i].len * h->n_smpl);
		l += b->gi[i].len * h->n_smpl;
	}
	return l;
}

int bcf_destroy(bcf1_t *b)
{
	int i;
	if (b == 0) return -1;
	free(b->str);
	for (i = 0; i < b->n_gi; ++i)
		free(b->gi[i].data);
	free(b->gi);
	free(b);
	return 0;
}

static inline void fmt_str(const char *p, kstring_t *s)
{
	if (*p == 0) kputc('.', s);
	else kputs(p, s);
}

void bcf_fmt_core(const bcf_hdr_t *h, bcf1_t *b, kstring_t *s)
{
	int i, j, x;
	s->l = 0;
	kputs(h->ns[b->tid], s); kputc('\t', s);
	kputw(b->pos + 1, s); kputc('\t', s);
	fmt_str(b->str, s); kputc('\t', s);
	fmt_str(b->ref, s); kputc('\t', s);
	fmt_str(b->alt, s); kputc('\t', s);
	kputw(b->qual, s); kputc('\t', s);
	fmt_str(b->flt, s); kputc('\t', s);
	fmt_str(b->info, s);
	if (b->fmt[0]) {
		kputc('\t', s);
		fmt_str(b->fmt, s);
	}
	x = b->n_alleles * (b->n_alleles + 1) / 2;
	for (j = 0; j < h->n_smpl; ++j) {
		kputc('\t', s);
		for (i = 0; i < b->n_gi; ++i) {
			if (i) kputc(':', s);
			if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
				uint8_t *d = (uint8_t*)b->gi[i].data + j * x;
				int k;
				for (k = 0; k < x; ++k) {
					if (k > 0) kputc(',', s);
					kputw(d[k], s);
				}
			} else if (b->gi[i].fmt == bcf_str2int("DP", 2)) {
				kputw(((uint16_t*)b->gi[i].data)[j], s);
			} else if (b->gi[i].fmt == bcf_str2int("GQ", 2)) {
				kputw(((uint8_t*)b->gi[i].data)[j], s);
			} else if (b->gi[i].fmt == bcf_str2int("GT", 2)) {
				int y = ((uint8_t*)b->gi[i].data)[j];
				kputc('0' + (y>>3&7), s);
				kputc("/|"[y>>6&1], s);
				kputc('0' + (y&7), s);
			}
		}
	}
}

char *bcf_fmt(const bcf_hdr_t *h, bcf1_t *b)
{
	kstring_t s;
	s.l = s.m = 0; s.s = 0;
	bcf_fmt_core(h, b, &s);
	return s.s;
}
