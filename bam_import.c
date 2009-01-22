#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include "bam.h"
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

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

struct __tamFile_t {
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	uint64_t n_lines;
};

char **bam_load_pos(const char *fn, int *_n)
{
	char **list = 0, *s;
	int n = 0, dret, m = 0, c;
	gzFile fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	kstream_t *ks;
	kstring_t *str;
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) > 0) {
		if (n == m) {
			m = m? m << 1 : 16;
			list = (char**)realloc(list, m * sizeof(char*));
		}
		s = list[n++] = (char*)calloc(str->l + 5, 1);
		strcpy(s, str->s);
		s += str->l + 1;
		ks_getuntil(ks, 0, str, &dret);
		*((uint32_t*)s) = atoi(str->s);
		if (dret != '\n')
			while ((c = ks_getc(fp)) >= 0 && c != '\n');
	}
	ks_destroy(ks);
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
	int c, dret, ret;
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	kh_ref_t *hash;
	khiter_t k;
	hash = kh_init(ref);
	fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	assert(fp);
	ks = ks_init(fp);
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		char *s = strdup(str->s);
		int len, i;
		i = kh_size(hash);
		ks_getuntil(ks, 0, str, &dret);
		len = atoi(str->s);
		k = kh_put(ref, hash, s, &ret);
		kh_value(hash, k) = (uint64_t)len<<32 | i;
		if (dret != '\n')
			while ((c = ks_getc(ks)) != '\n' && c != -1);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	fprintf(stderr, "[sam_header_read2] %d sequences loaded.\n", kh_size(hash));
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
	int x = header->l_text, y = header->l_text + str->l + 2; // 2 = 1 byte dret + 1 byte null
	kroundup32(x); kroundup32(y);
	if (x < y) header->text = (char*)realloc(header->text, y);
	strncpy(header->text + header->l_text, str->s, str->l+1); // we cannot use strcpy() here.
	header->l_text += str->l + 1;
	header->text[header->l_text] = 0;
}
int sam_read1(tamFile fp, bam_header_t *header, bam1_t *b)
{
	int ret, doff, doff0, dret;
	bam1_core_t *c = &b->core;
	kstring_t *str = fp->str;
	kstream_t *ks = fp->ks;

	while ((ret = ks_getuntil(fp->ks, 0, str, &dret)) >= 0 && str->s[0] == '@') { // skip header
		str->s[str->l] = dret; // note that str->s is NOT null terminated!!
		append_text(header, str);
		if (dret != '\n') {
			ret = ks_getuntil(fp->ks, '\n', str, &dret);
			str->s[str->l] = '\n'; // NOT null terminated!!
			append_text(header, str);
		}
		++fp->n_lines;
	}
	while (ret == 0) ret = ks_getuntil(fp->ks, 0, str, &dret); // special consideration for "\r\n"
	if (ret < 0) return -1;
	++fp->n_lines;
	doff = 0;

	{ // name
		c->l_qname = strlen(str->s) + 1;
		memcpy(alloc_data(b, doff + c->l_qname) + doff, str->s, c->l_qname);
		doff += c->l_qname;
	}
	{ // flag, tid, pos, qual
		ret = ks_getuntil(ks, 0, str, &dret); c->flag = atoi(str->s);
		ret = ks_getuntil(ks, 0, str, &dret); c->tid = bam_get_tid(header, str->s);
		ret = ks_getuntil(ks, 0, str, &dret); c->pos = isdigit(str->s[0])? atoi(str->s) - 1 : -1;
		ret = ks_getuntil(ks, 0, str, &dret); c->qual = isdigit(str->s[0])? atoi(str->s) : 0;
		if (ret < 0) return -2;
	}
	{ // cigar
		char *s, *t;
		int i, op;
		long x;
		c->n_cigar = 0;
		if (ks_getuntil(ks, 0, str, &dret) < 0) return -3;
		if (str->s[0] != '*') {
			for (s = str->s; *s; ++s) {
				if (isalpha(*s)) ++c->n_cigar;
				else if (!isdigit(*s)) parse_error(fp->n_lines, "invalid CIGAR character");
			}
			b->data = alloc_data(b, doff + c->n_cigar * 4);
			for (i = 0, s = str->s; i != c->n_cigar; ++i) {
				x = strtol(s, &t, 10);
				op = toupper(*t);
				if (op == 'M') op = BAM_CMATCH;
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
		}
	}
	{ // mtid, mpos, isize
		ret = ks_getuntil(ks, 0, str, &dret); c->mtid = strcmp(str->s, "=")? bam_get_tid(header, str->s) : c->tid;
		ret = ks_getuntil(ks, 0, str, &dret); c->mpos = isdigit(str->s[0])? atoi(str->s) - 1 : -1;
		ret = ks_getuntil(ks, 0, str, &dret); c->isize = (str->s[0] == '-' || isdigit(str->s[0]))? atoi(str->s) : 0;
		if (ret < 0) return -4;
	}
	{ // seq and qual
		int i;
		uint8_t *p;
		if (ks_getuntil(ks, 0, str, &dret) < 0) return -5; // seq
		c->l_qseq = strlen(str->s);
		if (c->n_cigar && c->l_qseq != (int32_t)bam_cigar2qlen(c, bam1_cigar(b)))
			parse_error(fp->n_lines, "CIGAR and sequence length are inconsistent");
		p = (uint8_t*)alloc_data(b, doff + c->l_qseq + (c->l_qseq+1)/2) + doff;
		bzero(p, (c->l_qseq+1)/2);
		for (i = 0; i < c->l_qseq; ++i)
			p[i/2] |= bam_nt16_table[(int)str->s[i]] << 4*(1-i%2);
		if (ks_getuntil(ks, 0, str, &dret) < 0) return -6; // qual
		if (c->l_qseq != strlen(str->s))
			parse_error(fp->n_lines, "sequence and quality are inconsistent");
		p += (c->l_qseq+1)/2;
		for (i = 0; i < c->l_qseq; ++i) p[i] = str->s[i] - 33;
		doff += c->l_qseq + (c->l_qseq+1)/2;
	}
	doff0 = doff;
	if (dret != '\n' && dret != '\r') { // aux
		while (ks_getuntil(ks, 0, str, &dret) >= 0) {
			uint8_t *s, type, key[2];
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
			} else parse_error(fp->n_lines, "unrecognized type");
			if (dret == '\n' || dret == '\r') break;
		}
	}
	b->l_aux = doff - doff0;
	b->data_len = doff;
	return 0;
}

tamFile sam_open(const char *fn)
{
	tamFile fp;
	fp = (tamFile)calloc(1, sizeof(struct __tamFile_t));
	fp->str = (kstring_t*)calloc(1, sizeof(kstring_t));
	fp->fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	fp->ks = ks_init(fp->fp);
	fp->n_lines = 0;
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

static void taf2baf_core(const char *fntaf, const char *fnbaf, bam_header_t *header)
{
	bamFile fpbaf;
	bam1_t *b;
	tamFile fp;
	int ret;

	b = (bam1_t*)calloc(1, sizeof(bam1_t));
	fpbaf = (strcmp(fnbaf, "-") == 0)? bam_dopen(fileno(stdout), "w") : bam_open(fnbaf, "w");
	fp = sam_open(fntaf);
	ret = sam_read1(fp, header, b);
	bam_header_write(fpbaf, header);
	if (ret >= 0) {
		bam_write1(fpbaf, b);
		while (sam_read1(fp, header, b) >= 0) bam_write1(fpbaf, b);
	}
	bam_close(fpbaf);
	free(b->data); free(b);
	sam_close(fp);
}

int bam_taf2baf(int argc, char *argv[])
{
	int c;
	bam_header_t *header;

	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 3 > argc) {
		fprintf(stderr, "Usage: bamtk import <in.ref_list> <in.sam> <out.bam>\n");
		return 1;
	}
	header = sam_header_read2(argv[optind]);
	taf2baf_core(argv[optind+1], argv[optind+2], header);
	bam_header_destroy(header);
	return 0;
}
