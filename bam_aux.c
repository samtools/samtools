#include <ctype.h>
#include "bam.h"
#include "khash.h"
KHASH_MAP_INIT_INT(aux, uint8_t*)
KHASH_MAP_INIT_STR(s, int)

void bam_init_header_hash(bam_header_t *header)
{
	if (header->hash == 0) {
		int ret, i;
		khiter_t iter;
		khash_t(s) *h;
		header->hash = h = kh_init(s);
		for (i = 0; i < header->n_targets; ++i) {
			iter = kh_put(s, h, header->target_name[i], &ret);
			kh_value(h, iter) = i;
		}
	}
}

void bam_destroy_header_hash(bam_header_t *header)
{
	if (header->hash)
		kh_destroy(s, (khash_t(s)*)header->hash);
}

int32_t bam_get_tid(const bam_header_t *header, const char *seq_name)
{
	khint_t k;
	khash_t(s) *h = (khash_t(s)*)header->hash;
	k = kh_get(s, h, seq_name);
	return k == kh_end(h)? -1 : kh_value(h, k);
}

void bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *begin, int *end)
{
	char *s, *p;
	int i, l, k;
	khiter_t iter;
	khash_t(s) *h;

	bam_init_header_hash(header);
	h = (khash_t(s)*)header->hash;
	
	l = strlen(str);
	p = s = (char*)malloc(l+1);
	/* squeeze out "," */
	for (i = k = 0; i != l; ++i)
		if (str[i] != ',' && !isspace(str[i])) s[k++] = str[i];
	s[k] = 0;
	for (i = 0; i != k; ++i) if (s[i] == ':') break;
	s[i] = 0;
	iter = kh_get(s, h, s); /* get the ref_id */
	if (iter == kh_end(h)) { // name not found
		*ref_id = -1; free(s);
		return;
	}
	*ref_id = kh_value(h, iter);
	if (i == k) { /* dump the whole sequence */
		*begin = 0; *end = 1<<29; free(s);
		return;
	}
	for (p = s + i + 1; i != k; ++i) if (s[i] == '-') break;
	*begin = atoi(p);
	if (i < k) {
		p = s + i + 1;
		*end = atoi(p);
	} else *end = 1<<29;
	if (*begin > 0) --*begin;
	assert(*begin <= *end);
	free(s);
}

void bam_aux_init(bam1_t *b)
{
	khash_t(aux) *h;
	uint8_t *s;
	if (b->hash == 0) {
		h = kh_init(aux);
		b->hash = h;
	} else {
		h = (khash_t(aux)*)b->hash;
		kh_clear(aux, h);
	}
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		uint32_t x = (uint32_t)s[0]<<8 | s[1];
		int ret, type;
		khint_t k;
		s += 2; type = toupper(*s); ++s;
		k = kh_put(aux, h, x, &ret);
		kh_value(h, k) = s;
		if (type == 'C') ++s;
		else if (type == 'S') s += 2;
		else if (type == 'I') s += 4;
		else if (type == 'F') s += 4;
		else if (type == 'Z') { while (*s) putchar(*s++); ++s; }
	}
}
void bam_aux_destroy(bam1_t *b)
{
	khash_t(aux) *h = (khash_t(aux)*)b->hash;
	kh_destroy(aux, h);
	b->hash = 0;
}
static uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2])
{
	uint32_t x = (uint32_t)tag[0]<<8 | tag[1];
	khint_t k;
	khash_t(aux) *h;
	if (b->hash == 0) bam_aux_init(b);
	h = (khash_t(aux)*)b->hash;
	k = kh_get(aux, h, x);
	if (k == kh_end(h)) return 0;
	return kh_value(h, k);
}
int32_t bam_aux_geti(bam1_t *b, const char tag[2], int *err)
{
	int type;
	uint8_t *s = bam_aux_get_core(b, tag);
	*err = 0;
	if (s == 0) { *err = -1; return 0; }
	type = *s++;
	if (type == 'c') return (int32_t)*(int8_t*)s;
	else if (type == 'C') return (int32_t)*(uint8_t*)s;
	else if (type == 's') return (int32_t)*(int16_t*)s;
	else if (type == 'S') return (int32_t)*(uint16_t*)s;
	else if (type == 'i' || type == 'I') return *(int32_t*)s;
	else { *err = -2; return 0; }
}
float bam_aux_getf(bam1_t *b, const char tag[2], int *err)
{
	int type;
	uint8_t *s = bam_aux_get_core(b, tag);
	*err = 0;
	type = *s++;
	if (s == 0) { *err = -1; return 0; }
	if (type == 'f') return *(float*)s;
	else { *err = -2; return 0; }
}
char bam_aux_getc(bam1_t *b, const char tag[2], int *err)
{
	int type;
	uint8_t *s = bam_aux_get_core(b, tag);
	*err = 0;
	type = *s++;
	if (s == 0) { *err = -1; return 0; }
	if (type == 'c') return *(char*)s;
	else { *err = -2; return 0; }
}
char *bam_aux_getZH(bam1_t *b, const char tag[2], int *err)
{
	int type;
	uint8_t *s = bam_aux_get_core(b, tag);
	*err = 0;
	type = *s++;
	if (s == 0) { *err = -1; return 0; }
	if (type == 'Z' || type == 'H') return (char*)s;
	else { *err = -2; return 0; }
}
