#include "bcf.h"
#include "kstring.h"
#include "khash.h"
KHASH_MAP_INIT_STR(str2id, int)

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
	int i, j, k, *z, n_smpl = b->n_smpl;
	if (b->n_alleles <= n) return -1;
	if (n > 1) {
		for (p = b->alt, k = 1; *p; ++p)
			if (*p == ',' && ++k == n) break;
		*p = '\0';
	} else p = b->alt, *p = '\0';
	++p;
	memmove(p, b->flt, b->str + b->l_str - b->flt);
	b->l_str -= b->flt - p;
	z = alloca(sizeof(int) / 2 * n * (n+1));
	for (i = k = 0; i < n; ++i)
		for (j = 0; j < n - i; ++j)
			z[k++] = i * b->n_alleles + j;
	for (i = 0; i < b->n_gi; ++i) {
		bcf_ginfo_t *g = b->gi + i;
		if (g->fmt == bcf_str2int("PL", 2)) {
			int l, x = b->n_alleles * (b->n_alleles + 1) / 2;
			uint8_t *d = (uint8_t*)g->data;
			g->len = n * (n + 1) / 2;
			for (l = k = 0; l < n_smpl; ++l) {
				uint8_t *dl = d + l * x;
				for (j = 0; j < g->len; ++j) d[k++] = dl[z[j]];
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
