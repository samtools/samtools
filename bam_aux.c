#include <ctype.h>
#include "bam.h"
#include "khash.h"
KHASH_MAP_INIT_STR(s, int)

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

uint8_t *bam_aux_get(bam1_t *b, const char tag[2])
{
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		int type, x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return s;
		type = toupper(*s); ++s;
		if (type == 'C') ++s;
		else if (type == 'S') s += 2;
		else if (type == 'I' || type == 'F') s += 4;
		else if (type == 'D') s += 8;
		else if (type == 'Z' || type == 'H') { while (*s) putchar(*s++); ++s; }
	}
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

char bam_aux_getCSi(bam1_t *b, int i)
{
	uint8_t *c = bam_aux_get(b, "CS");
	char *cs = NULL;

	// return the base if the tag was not found
	if(0 == c) return 0;

	cs = bam_aux2Z(c);
	// adjust for strandedness and leading adaptor
	if(bam1_strand(b)) i = strlen(cs) - 1 - i;
	else i++;
	return cs[i];
}

char bam_aux_getCQi(bam1_t *b, int i)
{
	uint8_t *c = bam_aux_get(b, "CQ");
	char *cq = NULL;
	
	// return the base if the tag was not found
	if(0 == c) return 0;

	cq = bam_aux2Z(c);
	// adjust for strandedness
	if(bam1_strand(b)) i = strlen(cq) - 1 - i;
	return cq[i];
}

char bam_aux_nt2int(char a)
{
	switch(toupper(a)) {
		case 'A':
			return 0;
			break;
		case 'C':
			return 1;
			break;
		case 'G':
			return 2;
			break;
		case 'T':
			return 3;
			break;
		default:
			return 4;
			break;
	}
}

char bam_aux_ntnt2cs(char a, char b)
{
	a = bam_aux_nt2int(a);
	b = bam_aux_nt2int(b);
	if(4 == a || 4 == b) return '4';
	return "0123"[(int)(a ^ b)];
}

char bam_aux_getCEi(bam1_t *b, int i)
{
	int cs_i;
	uint8_t *c = bam_aux_get(b, "CS");
	char *cs = NULL;
	char prev_b, cur_b;
	char cur_color, cor_color;

	// return the base if the tag was not found
	if(0 == c) return 0;
	
	cs = bam_aux2Z(c);

	// adjust for strandedness and leading adaptor
	if(bam1_strand(b)) { //reverse strand
		cs_i = strlen(cs) - 1 - i;
		// get current color
		cur_color = cs[cs_i];
		// get previous base
		prev_b = (0 == cs_i) ? cs[0] : bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i+1)];
		// get current base
		cur_b = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)]; 
	}
	else {
		cs_i=i+1;
		// get current color
		cur_color = cs[cs_i];
		// get previous base
		prev_b = (0 == i) ? cs[0] : bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i-1)];
		// get current base
		cur_b = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)];
	}

	// corrected color
	cor_color = bam_aux_ntnt2cs(prev_b, cur_b);

	if(cur_color == cor_color) { 
		return '-';
	}
	else {
		return cur_color;
	}
}
