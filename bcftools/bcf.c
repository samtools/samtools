#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "kstring.h"
#include "bcf.h"

int bcf_hdr_read(bcf_t *b);

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
		b->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_fdopen(fileno(stdout), "w");
	} else {
		b->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_fdopen(fileno(stdin), "r");
	}
	b->fp->owned_file = 1;
	if (strchr(mode, 'r')) bcf_hdr_read(b);
	return b;
}

int bcf_close(bcf_t *b)
{
	int ret;
	if (b == 0) return 0;
	ret = bgzf_close(b->fp);
	free(b->h.name); free(b->h.sname); free(b->h.txt); free(b->h.ns); free(b->h.sns);
	free(b);
	return ret;
}

int bcf_hdr_write(bcf_t *b)
{
	bcf_hdr_t *h;
	if (b == 0) return -1;
	h = &b->h;
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

int bcf_hdr_cpy(bcf_hdr_t *h, const bcf_hdr_t *h0)
{
	*h = *h0;
	h->name = malloc(h->l_nm);
	h->sname = malloc(h->l_smpl);
	h->txt = malloc(h->l_txt);
	memcpy(h->name, h0->name, h->l_nm);
	memcpy(h->sname, h0->sname, h->l_smpl);
	memcpy(h->txt, h0->txt, h->l_txt);
	bcf_hdr_sync(h);
	return 0;
}

int bcf_hdr_read(bcf_t *b)
{
	uint8_t magic[4];
	bcf_hdr_t *h;
	if (b == 0) return -1;
	bcf_hdr_clear(&b->h);
	h = &b->h;
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
	bcf_hdr_sync(&b->h);
	return 16 + h->l_nm + h->l_smpl + h->l_txt;
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
	b->ns = cnt_null(b->l_nm, b->name, &b->n_ref);
	b->sns = cnt_null(b->l_smpl, b->sname, &b->n_smpl);
	return 0;
}

#define char2int(s) (((int)s[0])<<8|s[1])

int bcf_sync(int n_smpl, bcf1_t *b)
{
	char *p, *tmp[5], *s;
	int i, n;
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
	for (p = s = b->fmt, n = 0; *p; ++p) {
		if (*p == ':' || *(p+1) == 0) {
			char *q = *p == ':'? p : p + 1;
			if ((q - s) != 2) return -2;
			b->gi[n].fmt = char2int(s);
			s = q;
		}
	}
	// set gi[i].len
	for (i = 0; i < b->n_gi; ++i) {
		if (b->gi[i].fmt == char2int("PL")) {
			b->gi[i].len = b->n_alleles * (b->n_alleles + 1) / 2;
		} else if (b->gi[i].fmt == char2int("DP") || b->gi[i].fmt == char2int("HQ")) {
			b->gi[i].len = 2;
		} else if (b->gi[i].fmt == char2int("GQ") || b->gi[i].fmt == char2int("GT")) {
			b->gi[i].len = 1;
		} else if (b->gi[i].fmt == char2int("GL")) {
			b->gi[i].len = 4;
		}
		b->gi[i].data = realloc(b->gi[i].data, n_smpl * b->gi[i].len);
	}
	return 0;
}

int bcf_write(bcf_t *bp, const bcf1_t *b)
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
		bgzf_write(bp->fp, b->gi[i].data, b->gi[i].len * bp->h.n_smpl);
		l += b->gi[i].len * bp->h.n_smpl;
	}
	return l;
}

int bcf_read(bcf_t *bp, bcf1_t *b)
{
	int i, l = 0;
	uint32_t x;
	if (b == 0) return -1;
	if (bgzf_read(bp->fp, &b->tid, 4) == 0) return 0;
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
	bcf_sync(bp->h.n_smpl, b);
	for (i = 0; i < b->n_gi; ++i) {
		bgzf_read(bp->fp, b->gi[i].data, b->gi[i].len * bp->h.n_smpl);
		l += b->gi[i].len * bp->h.n_smpl;
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

char *bcf_fmt(bcf_t *bp, bcf1_t *b)
{
	kstring_t s;
	int i, j, x;
	memset(&s, 0, sizeof(kstring_t));
	kputs(bp->h.ns[b->tid], &s); kputc('\t', &s);
	kputw(b->pos + 1, &s); kputc('\t', &s);
	fmt_str(b->str, &s); kputc('\t', &s);
	fmt_str(b->ref, &s); kputc('\t', &s);
	fmt_str(b->alt, &s); kputc('\t', &s);
	kputw(b->qual, &s); kputc('\t', &s);
	fmt_str(b->flt, &s); kputc('\t', &s);
	fmt_str(b->info, &s);
	if (b->fmt[0]) {
		kputc('\t', &s);
		fmt_str(b->fmt, &s);
	}
	x = b->n_alleles * (b->n_alleles + 1) / 2;
	for (j = 0; j < bp->h.n_smpl; ++j) {
		kputc('\t', &s);
		for (i = 0; i < b->n_gi; ++i) {
			if (i) kputc(':', &s);
			if (b->gi[i].fmt == char2int("PL")) {
				uint8_t *d = (uint8_t*)b->gi[i].data + j * x;
				int k;
				for (k = 0; k < x; ++k) {
					if (k > 0) kputc(',', &s);
					kputw(d[k], &s);
				}
			} else if (b->gi[i].fmt == char2int("DP")) {
				kputw(((uint16_t*)b->gi[i].data)[j], &s);
			} else if (b->gi[i].fmt == char2int("GQ")) {
				kputw(((uint8_t*)b->gi[i].data)[j], &s);
			}
		}
	}
	return s.s;
}
