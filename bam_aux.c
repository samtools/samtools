#include <ctype.h>
#include "bam.h"
#include "khash.h"
typedef char *str_p;
KHASH_MAP_INIT_STR(s, int)
KHASH_MAP_INIT_STR(r2l, str_p)

void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data)
{
	int ori_len = b->data_len;
	b->data_len += 3 + len;
	b->l_aux += 3 + len;
	if (b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	b->data[ori_len] = tag[0]; b->data[ori_len + 1] = tag[1];
	b->data[ori_len + 2] = type;
	memcpy(b->data + ori_len + 3, data, len);
}

uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2])
{
	return bam_aux_get(b, tag);
}

#define __skip_tag(s) do { \
		int type = toupper(*(s));										\
		++(s);															\
		if (type == 'C' || type == 'A') ++(s);							\
		else if (type == 'S') (s) += 2;									\
		else if (type == 'I' || type == 'F') (s) += 4;					\
		else if (type == 'D') (s) += 8;									\
		else if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
	} while (0)

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return s;
		__skip_tag(s);
	}
	return 0;
}
// s MUST BE returned by bam_aux_get()
int bam_aux_del(bam1_t *b, uint8_t *s)
{
	uint8_t *p, *aux;
	aux = bam1_aux(b);
	p = s - 2;
	__skip_tag(s);
	memmove(p, s, b->l_aux - (s - aux));
	b->data_len -= s - p;
	b->l_aux -= s - p;
	return 0;
}

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

int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *begin, int *end)
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
		return -1;
	}
	*ref_id = kh_value(h, iter);
	if (i == k) { /* dump the whole sequence */
		*begin = 0; *end = 1<<29; free(s);
		return -1;
	}
	for (p = s + i + 1; i != k; ++i) if (s[i] == '-') break;
	*begin = atoi(p);
	if (i < k) {
		p = s + i + 1;
		*end = atoi(p);
	} else *end = 1<<29;
	if (*begin > 0) --*begin;
	free(s);
	if (*begin > *end) {
		fprintf(stderr, "[bam_parse_region] invalid region.\n");
		return -1;
	}
	return 0;
}

int32_t bam_aux2i(const uint8_t *s)
{
	int type;
	if (s == 0) return 0;
	type = *s++;
	if (type == 'c') return (int32_t)*(int8_t*)s;
	else if (type == 'C') return (int32_t)*(uint8_t*)s;
	else if (type == 's') return (int32_t)*(int16_t*)s;
	else if (type == 'S') return (int32_t)*(uint16_t*)s;
	else if (type == 'i' || type == 'I') return *(int32_t*)s;
	else return 0;
}

float bam_aux2f(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0.0;
	if (type == 'f') return *(float*)s;
	else return 0.0;
}

double bam_aux2d(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0.0;
	if (type == 'd') return *(double*)s;
	else return 0.0;
}

char bam_aux2A(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0;
	if (type == 'A') return *(char*)s;
	else return 0;
}

char *bam_aux2Z(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0;
	if (type == 'Z' || type == 'H') return (char*)s;
	else return 0;
}
