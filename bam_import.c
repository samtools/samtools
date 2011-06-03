#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#ifdef _WIN32
#include <fcntl.h>
#endif
#include "kstring.h"
#include "bam.h"
#include "sam_header.h"
#include "kseq.h"
#include "khash.h"

KSTREAM_INIT(gzFile, gzread, 8192)
KHASH_MAP_INIT_STR(ref, uint64_t)

void bam_init_header_hash(bam_header_t *header);
void bam_destroy_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

unsigned short bam_char2flag_table[256] = {
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,BAM_FREAD1,BAM_FREAD2,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	BAM_FPROPER_PAIR,0,BAM_FMREVERSE,0, 0,BAM_FMUNMAP,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, BAM_FDUP,0,BAM_FQCFAIL,0, 0,0,0,0, 0,0,0,0,
	BAM_FPAIRED,0,BAM_FREVERSE,BAM_FSECONDARY, 0,BAM_FUNMAP,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
};

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

struct __tamFile_t {
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	uint64_t n_lines;
	int is_first;
};

char **__bam_get_lines(const char *fn, int *_n) // for bam_plcmd.c only
{
	char **list = 0, *s;
	int n = 0, dret, m = 0;
	gzFile fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	kstream_t *ks;
	kstring_t *str;
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	ks = ks_init(fp);
	while (ks_getuntil(ks, '\n', str, &dret) > 0) {
		if (n == m) {
			m = m? m << 1 : 16;
			list = (char**)realloc(list, m * sizeof(char*));
		}
		if (str->s[str->l-1] == '\r')
			str->s[--str->l] = '\0';
		s = list[n++] = (char*)calloc(str->l + 1, 1);
		strcpy(s, str->s);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	*_n = n;
	return list;
}

static bam_header_t *hash2header(const kh_ref_t *hash)
{
	bam_header_t *header;
	khiter_t k;
	header = bam_header_init();
	header->n_targets = kh_size(hash);
	header->target_name = (char**)calloc(kh_size(hash), sizeof(char*));
	header->target_len = (uint32_t*)calloc(kh_size(hash), 4);
	for (k = kh_begin(hash); k != kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			int i = (int)kh_value(hash, k);
			header->target_name[i] = (char*)kh_key(hash, k);
			header->target_len[i] = kh_value(hash, k)>>32;
		}
	}
	bam_init_header_hash(header);
	return header;
}
bam_header_t *sam_header_read2(const char *fn)
{
	bam_header_t *header;
	int c, dret, ret, error = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	kh_ref_t *hash;
	khiter_t k;
	if (fn == 0) return 0;
	fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	if (fp == 0) return 0;
	hash = kh_init(ref);
	ks = ks_init(fp);
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	while (ks_getuntil(ks, 0, str, &dret) > 0) {
		char *s = strdup(str->s);
		int len, i;
		i = kh_size(hash);
		ks_getuntil(ks, 0, str, &dret);
		len = atoi(str->s);
		k = kh_put(ref, hash, s, &ret);
		if (ret == 0) {
			fprintf(stderr, "[sam_header_read2] duplicated sequence name: %s\n", s);
			error = 1;
		}
		kh_value(hash, k) = (uint64_t)len<<32 | i;
		if (dret != '\n')
			while ((c = ks_getc(ks)) != '\n' && c != -1);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	fprintf(stderr, "[sam_header_read2] %d sequences loaded.\n", kh_size(hash));
	if (error) return 0;
	header = hash2header(hash);
	kh_destroy(ref, hash);
	return header;
}
static inline uint8_t *alloc_data(bam1_t *b, int size)
{
	if (b->m_data < size) {
		b->m_data = size;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	return b->data;
}
static inline void parse_error(int64_t n_lines, const char * __restrict msg)
{
	fprintf(stderr, "Parse error at line %lld: %s\n", (long long)n_lines, msg);
	abort();
}
static inline void append_text(bam_header_t *header, kstring_t *str)
{
	size_t x = header->l_text, y = header->l_text + str->l + 2; // 2 = 1 byte dret + 1 byte null
	kroundup32(x); kroundup32(y);
	if (x < y) 
    {
        header->n_text = y;
        header->text = (char*)realloc(header->text, y);
        if ( !header->text ) 
        {
            fprintf(stderr,"realloc failed to alloc %ld bytes\n", y);
            abort();
        }
    }
    // Sanity check
    if ( header->l_text+str->l+1 >= header->n_text )
    {
        fprintf(stderr,"append_text FIXME: %ld>=%ld, x=%ld,y=%ld\n",  header->l_text+str->l+1,header->n_text,x,y);
        abort();
    }
	strncpy(header->text + header->l_text, str->s, str->l+1); // we cannot use strcpy() here.
	header->l_text += str->l + 1;
	header->text[header->l_text] = 0;
}

int sam_header_parse(bam_header_t *h)
{
	char **tmp;
	int i;
	free(h->target_len); free(h->target_name);
	h->n_targets = 0; h->target_len = 0; h->target_name = 0;
	if (h->l_text < 3) return 0;
	if (h->dict == 0) h->dict = sam_header_parse2(h->text);
	tmp = sam_header2list(h->dict, "SQ", "SN", &h->n_targets);
	if (h->n_targets == 0) return 0;
	h->target_name = calloc(h->n_targets, sizeof(void*));
	for (i = 0; i < h->n_targets; ++i)
		h->target_name[i] = strdup(tmp[i]);
	free(tmp);
	tmp = sam_header2list(h->dict, "SQ", "LN", &h->n_targets);
	h->target_len = calloc(h->n_targets, 4);
	for (i = 0; i < h->n_targets; ++i)
		h->target_len[i] = atoi(tmp[i]);
	free(tmp);
	return h->n_targets;
}

bam_header_t *sam_header_read(tamFile fp)
{
	int ret, dret;
	bam_header_t *header = bam_header_init();
	kstring_t *str = fp->str;
	while ((ret = ks_getuntil(fp->ks, KS_SEP_TAB, str, &dret)) >= 0 && str->s[0] == '@') { // skip header
		str->s[str->l] = dret; // note that str->s is NOT null terminated!!
		append_text(header, str);
		if (dret != '\n') {
			ret = ks_getuntil(fp->ks, '\n', str, &dret);
			str->s[str->l] = '\n'; // NOT null terminated!!
			append_text(header, str);
		}
		++fp->n_lines;
	}
	sam_header_parse(header);
	bam_init_header_hash(header);
	fp->is_first = 1;
	return header;
}

int sam_read1(tamFile fp, bam_header_t *header, bam1_t *b)
{
	int ret, doff, doff0, dret, z = 0;
	bam1_core_t *c = &b->core;
	kstring_t *str = fp->str;
	kstream_t *ks = fp->ks;

	if (fp->is_first) {
		fp->is_first = 0;
		ret = str->l;
	} else {
		do { // special consideration for empty lines
			ret = ks_getuntil(fp->ks, KS_SEP_TAB, str, &dret);
			if (ret >= 0) z += str->l + 1;
		} while (ret == 0);
	}
	if (ret < 0) return -1;
	++fp->n_lines;
	doff = 0;

	{ // name
		c->l_qname = strlen(str->s) + 1;
		memcpy(alloc_data(b, doff + c->l_qname) + doff, str->s, c->l_qname);
		doff += c->l_qname;
	}
	{ // flag
		long flag;
		char *s;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		flag = strtol((char*)str->s, &s, 0);
		if (*s) { // not the end of the string
			flag = 0;
			for (s = str->s; *s; ++s)
				flag |= bam_char2flag_table[(int)*s];
		}
		c->flag = flag;
	}
	{ // tid, pos, qual
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1; c->tid = bam_get_tid(header, str->s);
		if (c->tid < 0 && strcmp(str->s, "*")) {
			if (header->n_targets == 0) {
				fprintf(stderr, "[sam_read1] missing header? Abort!\n");
				exit(1);
			} else fprintf(stderr, "[sam_read1] reference '%s' is recognized as '*'.\n", str->s);
		}
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1; c->pos = isdigit(str->s[0])? atoi(str->s) - 1 : -1;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1; c->qual = isdigit(str->s[0])? atoi(str->s) : 0;
		if (ret < 0) return -2;
	}
	{ // cigar
		char *s, *t;
		int i, op;
		long x;
		c->n_cigar = 0;
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -3;
		z += str->l + 1;
		if (str->s[0] != '*') {
			for (s = str->s; *s; ++s) {
				if (isalpha(*s)) ++c->n_cigar;
				else if (!isdigit(*s)) parse_error(fp->n_lines, "invalid CIGAR character");
			}
			b->data = alloc_data(b, doff + c->n_cigar * 4);
			for (i = 0, s = str->s; i != c->n_cigar; ++i) {
				x = strtol(s, &t, 10);
				op = toupper(*t);
				if (op == 'M' || op == '=' || op == 'X') op = BAM_CMATCH;
				else if (op == 'I') op = BAM_CINS;
				else if (op == 'D') op = BAM_CDEL;
				else if (op == 'N') op = BAM_CREF_SKIP;
				else if (op == 'S') op = BAM_CSOFT_CLIP;
				else if (op == 'H') op = BAM_CHARD_CLIP;
				else if (op == 'P') op = BAM_CPAD;
				else parse_error(fp->n_lines, "invalid CIGAR operation");
				s = t + 1;
				bam1_cigar(b)[i] = x << BAM_CIGAR_SHIFT | op;
			}
			if (*s) parse_error(fp->n_lines, "unmatched CIGAR operation");
			c->bin = bam_reg2bin(c->pos, bam_calend(c, bam1_cigar(b)));
			doff += c->n_cigar * 4;
		} else {
			if (!(c->flag&BAM_FUNMAP)) {
				fprintf(stderr, "Parse warning at line %lld: mapped sequence without CIGAR\n", (long long)fp->n_lines);
				c->flag |= BAM_FUNMAP;
			}
			c->bin = bam_reg2bin(c->pos, c->pos + 1);
		}
	}
	{ // mtid, mpos, isize
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		c->mtid = strcmp(str->s, "=")? bam_get_tid(header, str->s) : c->tid;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		c->mpos = isdigit(str->s[0])? atoi(str->s) - 1 : -1;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		c->isize = (str->s[0] == '-' || isdigit(str->s[0]))? atoi(str->s) : 0;
		if (ret < 0) return -4;
	}
	{ // seq and qual
		int i;
		uint8_t *p = 0;
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -5; // seq
		z += str->l + 1;
		if (strcmp(str->s, "*")) {
			c->l_qseq = strlen(str->s);
			if (c->n_cigar && c->l_qseq != (int32_t)bam_cigar2qlen(c, bam1_cigar(b)))
				parse_error(fp->n_lines, "CIGAR and sequence length are inconsistent");
			p = (uint8_t*)alloc_data(b, doff + c->l_qseq + (c->l_qseq+1)/2) + doff;
			memset(p, 0, (c->l_qseq+1)/2);
			for (i = 0; i < c->l_qseq; ++i)
				p[i/2] |= bam_nt16_table[(int)str->s[i]] << 4*(1-i%2);
		} else c->l_qseq = 0;
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -6; // qual
		z += str->l + 1;
		if (strcmp(str->s, "*") && c->l_qseq != strlen(str->s))
			parse_error(fp->n_lines, "sequence and quality are inconsistent");
		p += (c->l_qseq+1)/2;
		if (strcmp(str->s, "*") == 0) for (i = 0; i < c->l_qseq; ++i) p[i] = 0xff;
		else for (i = 0; i < c->l_qseq; ++i) p[i] = str->s[i] - 33;
		doff += c->l_qseq + (c->l_qseq+1)/2;
	}
	doff0 = doff;
	if (dret != '\n' && dret != '\r') { // aux
		while (ks_getuntil(ks, KS_SEP_TAB, str, &dret) >= 0) {
			uint8_t *s, type, key[2];
			z += str->l + 1;
			if (str->l < 6 || str->s[2] != ':' || str->s[4] != ':')
				parse_error(fp->n_lines, "missing colon in auxiliary data");
			key[0] = str->s[0]; key[1] = str->s[1];
			type = str->s[3];
			s = alloc_data(b, doff + 3) + doff;
			s[0] = key[0]; s[1] = key[1]; s += 2; doff += 2;
			if (type == 'A' || type == 'a' || type == 'c' || type == 'C') { // c and C for backward compatibility
				s = alloc_data(b, doff + 2) + doff;
				*s++ = 'A'; *s = str->s[5];
				doff += 2;
			} else if (type == 'I' || type == 'i') {
				long long x;
				s = alloc_data(b, doff + 5) + doff;
				x = (long long)atoll(str->s + 5);
				if (x < 0) {
					if (x >= -127) {
						*s++ = 'c'; *(int8_t*)s = (int8_t)x;
						s += 1; doff += 2;
					} else if (x >= -32767) {
						*s++ = 's'; *(int16_t*)s = (int16_t)x;
						s += 2; doff += 3;
					} else {
						*s++ = 'i'; *(int32_t*)s = (int32_t)x;
						s += 4; doff += 5;
						if (x < -2147483648ll)
							fprintf(stderr, "Parse warning at line %lld: integer %lld is out of range.",
									(long long)fp->n_lines, x);
					}
				} else {
					if (x <= 255) {
						*s++ = 'C'; *s++ = (uint8_t)x;
						doff += 2;
					} else if (x <= 65535) {
						*s++ = 'S'; *(uint16_t*)s = (uint16_t)x;
						s += 2; doff += 3;
					} else {
						*s++ = 'I'; *(uint32_t*)s = (uint32_t)x;
						s += 4; doff += 5;
						if (x > 4294967295ll)
							fprintf(stderr, "Parse warning at line %lld: integer %lld is out of range.",
									(long long)fp->n_lines, x);
					}
				}
			} else if (type == 'f') {
				s = alloc_data(b, doff + 5) + doff;
				*s++ = 'f';
				*(float*)s = (float)atof(str->s + 5);
				s += 4; doff += 5;
			} else if (type == 'd') {
				s = alloc_data(b, doff + 9) + doff;
				*s++ = 'd';
				*(float*)s = (float)atof(str->s + 9);
				s += 8; doff += 9;
			} else if (type == 'Z' || type == 'H') {
				int size = 1 + (str->l - 5) + 1;
				if (type == 'H') { // check whether the hex string is valid
					int i;
					if ((str->l - 5) % 2 == 1) parse_error(fp->n_lines, "length of the hex string not even");
					for (i = 0; i < str->l - 5; ++i) {
						int c = toupper(str->s[5 + i]);
						if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')))
							parse_error(fp->n_lines, "invalid hex character");
					}
				}
				s = alloc_data(b, doff + size) + doff;
				*s++ = type;
				memcpy(s, str->s + 5, str->l - 5);
				s[str->l - 5] = 0;
				doff += size;
			} else if (type == 'B') {
				int32_t n = 0, Bsize, k = 0, size;
				char *p;
				if (str->l < 8) parse_error(fp->n_lines, "too few values in aux type B");
				Bsize = bam_aux_type2size(str->s[5]); // the size of each element
				for (p = (char*)str->s + 6; *p; ++p) // count the number of elements in the array
					if (*p == ',') ++n;
				p = str->s + 7; // now p points to the first number in the array
				size = 6 + Bsize * n; // total number of bytes allocated to this tag
				s = alloc_data(b, doff + 6 * Bsize * n) + doff; // allocate memory
				*s++ = 'B'; *s++ = str->s[5];
				memcpy(s, &n, 4); s += 4; // write the number of elements
				if (str->s[5] == 'c')      while (p < str->s + str->l) ((int8_t*)s)[k++]   = (int8_t)strtol(p, &p, 0),   ++p;
				else if (str->s[5] == 'C') while (p < str->s + str->l) ((uint8_t*)s)[k++]  = (uint8_t)strtol(p, &p, 0),  ++p;
				else if (str->s[5] == 's') while (p < str->s + str->l) ((int16_t*)s)[k++]  = (int16_t)strtol(p, &p, 0),  ++p; // FIXME: avoid unaligned memory
				else if (str->s[5] == 'S') while (p < str->s + str->l) ((uint16_t*)s)[k++] = (uint16_t)strtol(p, &p, 0), ++p;
				else if (str->s[5] == 'i') while (p < str->s + str->l) ((int32_t*)s)[k++]  = (int32_t)strtol(p, &p, 0),  ++p;
				else if (str->s[5] == 'I') while (p < str->s + str->l) ((uint32_t*)s)[k++] = (uint32_t)strtol(p, &p, 0), ++p;
				else if (str->s[5] == 'f') while (p < str->s + str->l) ((float*)s)[k++]    = (float)strtod(p, &p),       ++p;
				else parse_error(fp->n_lines, "unrecognized array type");
				s += Bsize * n; doff += size;
			} else parse_error(fp->n_lines, "unrecognized type");
			if (dret == '\n' || dret == '\r') break;
		}
	}
	b->l_aux = doff - doff0;
	b->data_len = doff;
	return z;
}

tamFile sam_open(const char *fn)
{
	tamFile fp;
	gzFile gzfp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "rb") : gzopen(fn, "rb");
	if (gzfp == 0) return 0;
	fp = (tamFile)calloc(1, sizeof(struct __tamFile_t));
	fp->str = (kstring_t*)calloc(1, sizeof(kstring_t));
	fp->fp = gzfp;
	fp->ks = ks_init(fp->fp);
	return fp;
}

void sam_close(tamFile fp)
{
	if (fp) {
		ks_destroy(fp->ks);
		gzclose(fp->fp);
		free(fp->str->s); free(fp->str);
		free(fp);
	}
}
