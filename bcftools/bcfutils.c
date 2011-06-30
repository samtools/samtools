#include <string.h>
#include <math.h>
#include "bcf.h"
#include "kstring.h"
#include "khash.h"
KHASH_MAP_INIT_STR(str2id, int)

// FIXME: valgrind report a memory leak in this function. Probably it does not get deallocated...
void *bcf_build_refhash(bcf_hdr_t *h)
{
	khash_t(str2id) *hash;
	int i, ret;
	hash = kh_init(str2id);
	for (i = 0; i < h->n_ref; ++i) {
		khint_t k;
		k = kh_put(str2id, hash, h->ns[i], &ret); // FIXME: check ret
		kh_val(hash, k) = i;
	}
	return hash;
}

void *bcf_str2id_init()
{
	return kh_init(str2id);
}

void bcf_str2id_destroy(void *_hash)
{
	khash_t(str2id) *hash = (khash_t(str2id)*)_hash;
	if (hash) kh_destroy(str2id, hash); // Note that strings are not freed.
}

void bcf_str2id_thorough_destroy(void *_hash)
{
	khash_t(str2id) *hash = (khash_t(str2id)*)_hash;
	khint_t k;
	if (hash == 0) return;
	for (k = 0; k < kh_end(hash); ++k)
		if (kh_exist(hash, k)) free((char*)kh_key(hash, k));
	kh_destroy(str2id, hash);
}

int bcf_str2id(void *_hash, const char *str)
{
	khash_t(str2id) *hash = (khash_t(str2id)*)_hash;
	khint_t k;
	if (!hash) return -1;
	k = kh_get(str2id, hash, str);
	return k == kh_end(hash)? -1 : kh_val(hash, k);
}

int bcf_str2id_add(void *_hash, const char *str)
{
	khint_t k;
	int ret;
	khash_t(str2id) *hash = (khash_t(str2id)*)_hash;
	if (!hash) return -1;
	k = kh_put(str2id, hash, str, &ret);
	if (ret == 0) return kh_val(hash, k);
	kh_val(hash, k) = kh_size(hash) - 1;
	return kh_val(hash, k);
}

int bcf_shrink_alt(bcf1_t *b, int n)
{
	char *p;
	int i, j, k, n_smpl = b->n_smpl;
	if (b->n_alleles <= n) return -1;
	// update ALT
	if (n > 1) {
		for (p = b->alt, k = 1; *p; ++p)
			if (*p == ',' && ++k == n) break;
		*p = '\0';
	} else p = b->alt, *p = '\0';
	++p;
	memmove(p, b->flt, b->str + b->l_str - b->flt);
	b->l_str -= b->flt - p;
	// update PL
	for (i = 0; i < b->n_gi; ++i) {
		bcf_ginfo_t *g = b->gi + i;
		if (g->fmt == bcf_str2int("PL", 2)) {
			int l, x = b->n_alleles * (b->n_alleles + 1) / 2;
			uint8_t *d = (uint8_t*)g->data;
			g->len = n * (n + 1) / 2;
			for (l = k = 0; l < n_smpl; ++l) {
				uint8_t *dl = d + l * x;
				for (j = 0; j < g->len; ++j) d[k++] = dl[j];
			}
		} // FIXME: to add GL
	}
	b->n_alleles = n;
	bcf_sync(b);
	return 0;
}

int bcf_gl2pl(bcf1_t *b)
{
	char *p;
	int i, n_smpl = b->n_smpl;
	bcf_ginfo_t *g;
	float *d0;
	uint8_t *d1;
	if (strstr(b->fmt, "PL")) return -1;
	if ((p = strstr(b->fmt, "GL")) == 0) return -1;
	*p = 'P';
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("GL", 2))
			break;
	g = b->gi + i;
	g->fmt = bcf_str2int("PL", 2);
	g->len /= 4; // 4 == sizeof(float)
	d0 = (float*)g->data; d1 = (uint8_t*)g->data;
	for (i = 0; i < n_smpl * g->len; ++i) {
		int x = (int)(-10. * d0[i] + .499);
		if (x > 255) x = 255;
		if (x < 0) x = 0;
		d1[i] = x;
	}
	return 0;
}
/* FIXME: this function will fail given AB:GTX:GT. BCFtools never
 * produces such FMT, but others may do. */
int bcf_fix_gt(bcf1_t *b)
{
	char *s;
	int i;
	uint32_t tmp;
	bcf_ginfo_t gt;
	// check the presence of the GT FMT
	if ((s = strstr(b->fmt, ":GT")) == 0) return 0; // no GT or GT is already the first
	if (s[3] != '\0' && s[3] != ':') return 0; // :GTX in fact
	tmp = bcf_str2int("GT", 2);
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == tmp) break;
	if (i == b->n_gi) return 0; // no GT in b->gi; probably a bug...
	gt = b->gi[i];
	// move GT to the first
	for (; i > 0; --i) b->gi[i] = b->gi[i-1];
	b->gi[0] = gt;
	memmove(b->fmt + 3, b->fmt, s + 1 - b->fmt);
	b->fmt[0] = 'G'; b->fmt[1] = 'T'; b->fmt[2] = ':';
	return 0;
}

int bcf_fix_pl(bcf1_t *b)
{
	int i;
	uint32_t tmp;
	uint8_t *PL, *swap;
	bcf_ginfo_t *gi;
	// pinpoint PL
	tmp = bcf_str2int("PL", 2);
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == tmp) break;
	if (i == b->n_gi) return 0;
	// prepare
	gi = b->gi + i;
	PL = (uint8_t*)gi->data;
	swap = alloca(gi->len);
	// loop through individuals
	for (i = 0; i < b->n_smpl; ++i) {
		int k, l, x;
		uint8_t *PLi = PL + i * gi->len;
		memcpy(swap, PLi, gi->len);
		for (k = x = 0; k < b->n_alleles; ++k)
			for (l = k; l < b->n_alleles; ++l)
				PLi[l*(l+1)/2 + k] = swap[x++];
	}
	return 0;
}

int bcf_smpl_covered(const bcf1_t *b)
{
	int i, j, n = 0;
	uint32_t tmp;
	bcf_ginfo_t *gi;
	// pinpoint PL
	tmp = bcf_str2int("PL", 2);
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == tmp) break;
	if (i == b->n_gi) return 0;
	// count how many samples having PL!=[0..0]
	gi = b->gi + i;
	for (i = 0; i < b->n_smpl; ++i) {
		uint8_t *PLi = ((uint8_t*)gi->data) + i * gi->len;
		for (j = 0; j < gi->len; ++j)
			if (PLi[j]) break;
		if (j < gi->len) ++n;
	}
	return n;
}

static void *locate_field(const bcf1_t *b, const char *fmt, int l)
{
	int i;
	uint32_t tmp;
	tmp = bcf_str2int(fmt, l);
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == tmp) break;
	return i == b->n_gi? 0 : b->gi[i].data;
}

int bcf_anno_max(bcf1_t *b)
{
	int k, max_gq, max_sp, n_het;
	kstring_t str;
	uint8_t *gt, *gq;
	int32_t *sp;
	max_gq = max_sp = n_het = 0;
	gt = locate_field(b, "GT", 2);
	if (gt == 0) return -1;
	gq = locate_field(b, "GQ", 2);
	sp = locate_field(b, "SP", 2);
	if (sp)
		for (k = 0; k < b->n_smpl; ++k)
			if (gt[k]&0x3f)
				max_sp = max_sp > (int)sp[k]? max_sp : sp[k];
	if (gq)
		for (k = 0; k < b->n_smpl; ++k)
			if (gt[k]&0x3f)
				max_gq = max_gq > (int)gq[k]? max_gq : gq[k];
	for (k = 0; k < b->n_smpl; ++k) {
		int a1, a2;
		a1 = gt[k]&7; a2 = gt[k]>>3&7;
		if ((!a1 && a2) || (!a2 && a1)) { // a het
			if (gq == 0) ++n_het;
			else if (gq[k] >= 20) ++n_het;
		}
	}
	if (n_het) max_sp -= (int)(4.343 * log(n_het) + .499);
	if (max_sp < 0) max_sp = 0;
	memset(&str, 0, sizeof(kstring_t));
	if (*b->info) kputc(';', &str);
	ksprintf(&str, "MXSP=%d;MXGQ=%d", max_sp, max_gq);
	bcf_append_info(b, str.s, str.l);
	free(str.s);
	return 0;
}

// FIXME: only data are shuffled; the header is NOT
int bcf_shuffle(bcf1_t *b, int seed)
{
	int i, j, *a;
	if (seed > 0) srand48(seed);
	a = malloc(b->n_smpl * sizeof(int));
	for (i = 0; i < b->n_smpl; ++i) a[i] = i;
	for (i = b->n_smpl; i > 1; --i) {
		int tmp;
		j = (int)(drand48() * i);
		tmp = a[j]; a[j] = a[i-1]; a[i-1] = tmp;
	}
	for (j = 0; j < b->n_gi; ++j) {
		bcf_ginfo_t *gi = b->gi + j;
		uint8_t *swap, *data = (uint8_t*)gi->data;
		swap = malloc(gi->len * b->n_smpl);
		for (i = 0; i < b->n_smpl; ++i)
			memcpy(swap + gi->len * a[i], data + gi->len * i, gi->len);
		free(gi->data);
		gi->data = swap;
	}
	free(a);
	return 0;
}

bcf_hdr_t *bcf_hdr_subsam(const bcf_hdr_t *h0, int n, char *const* samples, int *list)
{
	int i, ret, j;
	khint_t k;
	bcf_hdr_t *h;
	khash_t(str2id) *hash;
	kstring_t s;
	s.l = s.m = 0; s.s = 0;
	hash = kh_init(str2id);
	for (i = 0; i < h0->n_smpl; ++i) {
		k = kh_put(str2id, hash, h0->sns[i], &ret);
		kh_val(hash, k) = i;
	}
	for (i = j = 0; i < n; ++i) {
		k = kh_get(str2id, hash, samples[i]);
		if (k != kh_end(hash)) {
			list[j++] = kh_val(hash, k);
			kputs(samples[i], &s); kputc('\0', &s);
		}
	}
	if (j < n) fprintf(stderr, "<%s> %d samples in the list but not in BCF.", __func__, n - j);
	kh_destroy(str2id, hash);
	h = calloc(1, sizeof(bcf_hdr_t));
	*h = *h0;
	h->ns = 0; h->sns = 0;
	h->name = malloc(h->l_nm); memcpy(h->name, h0->name, h->l_nm);
	h->txt = calloc(1, h->l_txt + 1); memcpy(h->txt, h0->txt, h->l_txt);
	h->l_smpl = s.l; h->sname = s.s;
	bcf_hdr_sync(h);
	return h;
}

int bcf_subsam(int n_smpl, int *list, bcf1_t *b)
{
	int i, j;
	for (j = 0; j < b->n_gi; ++j) {
		bcf_ginfo_t *gi = b->gi + j;
		uint8_t *swap;
		swap = malloc(gi->len * b->n_smpl);
		for (i = 0; i < n_smpl; ++i)
			memcpy(swap + i * gi->len, (uint8_t*)gi->data + list[i] * gi->len, gi->len);
		free(gi->data);
		gi->data = swap;
	}
	b->n_smpl = n_smpl;
	return 0;
}

static int8_t nt4_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4, -1, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4, -1, 4, 4, 4,  4, 4, 4, 4
};

int bcf_gl10(const bcf1_t *b, uint8_t *gl)
{
	int a[4], k, l, map[4], k1, j, i;
	const bcf_ginfo_t *PL;
	char *s;
	if (b->ref[1] != 0 || b->n_alleles > 4) return -1; // ref is not a single base or >4 alleles
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) break;
	if (i == b->n_gi) return -1; // no PL
	PL = b->gi + i;
	a[0] = nt4_table[(int)b->ref[0]];
	if (a[0] > 3 || a[0] < 0) return -1; // ref is not A/C/G/T
	a[1] = a[2] = a[3] = -2; // -1 has a special meaning
	if (b->alt[0] == 0) return -1; // no alternate allele
	map[0] = map[1] = map[2] = map[3] = -2;
	map[a[0]] = 0;
	for (k = 0, s = b->alt, k1 = -1; k < 3 && *s; ++k, s += 2) {
		if (s[1] != ',' && s[1] != 0) return -1; // ALT is not single base
		a[k+1] = nt4_table[(int)*s];
		if (a[k+1] >= 0) map[a[k+1]] = k+1;
		else k1 = k + 1;
		if (s[1] == 0) break; // the end of the ALT string
	}
	for (k = 0; k < 4; ++k)
		if (map[k] < 0) map[k] = k1;
	for (i = 0; i < b->n_smpl; ++i) {
		const uint8_t *p = PL->data + i * PL->len; // the PL for the i-th individual
		uint8_t *g = gl + 10 * i;
		for (j = 0; j < PL->len; ++j)
			if (p[j]) break;
		for (k = j = 0; k < 4; ++k) {
			for (l = k; l < 4; ++l) {
				int t, x = map[k], y = map[l];
				if (x > y) t = x, x = y, y = t; // make sure x is the smaller
				g[j++] = p[y * (y+1) / 2 + x];
			}
		}
	}
	return 0;
}
