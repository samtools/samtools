#include "bcf.h"
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

int bcf_str2id_put(void *_hash, const char *str, int id)
{
	return 0;
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

