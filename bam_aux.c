#include <ctype.h>
#include "bam.h"

static inline int bam_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f' || x == 'F') return 4;
	else return 0;
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += bam_aux_type2size(type); \
	} while(0)


int bam_aux_drop_other(bam1_t *b, uint8_t *s)
{
	if (s) {
		uint8_t *p, *aux;
		aux = bam1_aux(b);
		p = s - 2;
		__skip_tag(s);
		memmove(aux, p, s - p);
		b->data_len -= bam_get_l_aux(b) - (s - p);
	} else {
		b->data_len -= bam_get_l_aux(b);
	}
	return 0;
}

int bam_parse_region(bam_hdr_t *header, const char *str, int *ref_id, int *beg, int *end)
{
	const char *name_lim = hts_parse_reg(str, beg, end);
	char *name = malloc(name_lim - str + 1);
	memcpy(name, str, name_lim - str);
	name[name_lim - str] = '\0';
	*ref_id = bam_name2id(header, name);
	free(name);
	if (*ref_id == -1) return -1;
	return *beg <= *end? 0 : -1;
}
