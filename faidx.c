#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "faidx.h"
#include "khash.h"

typedef struct {
	uint64_t len:32, line_len:16, line_blen:16;
	uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

#ifndef _NO_RAZF
#include "razf.h"
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#define RAZF FILE
#define razf_read(fp, buf, size) fread(buf, 1, size, fp)
#define razf_open(fn, mode) fopen(fn, mode)
#define razf_close(fp) fclose(fp)
#define razf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define razf_tell(fp) ftello(fp)
#endif

struct __faidx_t {
	RAZF *rz;
	int n, m;
	char **name;
	khash_t(s) *hash;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline void fai_insert_index(faidx_t *idx, const char *name, int len, int line_len, int line_blen, uint64_t offset)
{
	khint_t k;
	int ret;
	faidx1_t t;
	if (idx->n == idx->m) {
		idx->m = idx->m? idx->m<<1 : 16;
		idx->name = (char**)realloc(idx->name, sizeof(void*) * idx->m);
	}
	idx->name[idx->n] = strdup(name);
	k = kh_put(s, idx->hash, idx->name[idx->n], &ret);
	t.len = len; t.line_len = line_len; t.line_blen = line_blen; t.offset = offset;
	kh_value(idx->hash, k) = t;
	++idx->n;
}

faidx_t *fai_build_core(RAZF *rz)
{
	char c, *name;
	int l_name, m_name, ret;
	int len, line_len, line_blen, state;
	int l1, l2;
	faidx_t *idx;
	uint64_t offset;

	idx = (faidx_t*)calloc(1, sizeof(faidx_t));
	idx->hash = kh_init(s);
	name = 0; l_name = m_name = 0;
	len = line_len = line_blen = -1; state = 0; l1 = l2 = -1; offset = 0;
	while (razf_read(rz, &c, 1)) {
		if (c == '>') { // fasta header
			if (len >= 0)
				fai_insert_index(idx, name, len, line_len, line_blen, offset);
			l_name = 0;
			while ((ret = razf_read(rz, &c, 1)) != 0 && !isspace(c)) {
				if (m_name < l_name + 2) {
					m_name = l_name + 2;
					kroundup32(m_name);
					name = (char*)realloc(name, m_name);
				}
				name[l_name++] = c;
			}
			name[l_name] = '\0';
			assert(ret);
			if (c != '\n') while (razf_read(rz, &c, 1) && c != '\n');
			state = 1; len = 0;
			offset = razf_tell(rz);
		} else {
			if (state == 3) {
				fprintf(stderr, "[fai_build_core] inlined empty line is not allowed in sequence '%s'. Abort!\n", name);
				exit(1);
			}
			if (state == 2) state = 3;
			l1 = l2 = 0;
			do {
				++l1;
				if (isgraph(c)) ++l2;
			} while ((ret = razf_read(rz, &c, 1)) && c != '\n');
			if (state == 3 && l2) {
				fprintf(stderr, "[fai_build_core] different line length in sequence '%s'. Abort!\n", name);
				exit(1);
			}
			++l1; len += l2;
			if (l2 >= 0x10000) {
				fprintf(stderr, "[fai_build_core] line length exceeds 65535 in sequence '%s'. Abort!\n", name);
				exit(1);
			}
			if (state == 1) line_len = l1, line_blen = l2, state = 0;
			else if (state == 0) {
				if (l1 != line_len || l2 != line_blen) state = 2;
			}
		}
	}
	fai_insert_index(idx, name, len, line_len, line_blen, offset);
	free(name);
	return idx;
}

void fai_save(const faidx_t *fai, FILE *fp)
{
	khint_t k;
	int i;
	for (i = 0; i < fai->n; ++i) {
		faidx1_t x;
		k = kh_get(s, fai->hash, fai->name[i]);
		x = kh_value(fai->hash, k);
		fprintf(fp, "%s\t%d\t%lld\t%d\t%d\n", fai->name[i], (int)x.len, (long long)x.offset, (int)x.line_blen, (int)x.line_len);
	}
}

faidx_t *fai_read(FILE *fp)
{
	faidx_t *fai;
	char *buf, *p;
	int len, line_len, line_blen;
	long long offset;
	fai = (faidx_t*)calloc(1, sizeof(faidx_t));
	fai->hash = kh_init(s);
	buf = (char*)calloc(0x10000, 1);
	while (!feof(fp) && fgets(buf, 0x10000, fp)) {
		for (p = buf; *p && isgraph(*p); ++p);
		*p = 0; ++p;
		sscanf(p, "%d%lld%d%d", &len, &offset, &line_blen, &line_len);
		fai_insert_index(fai, buf, len, line_len, line_blen, offset);
	}
	free(buf);
	return fai;
}

void fai_destroy(faidx_t *fai)
{
	int i;
	for (i = 0; i < fai->n; ++i) free(fai->name[i]);
	free(fai->name);
	kh_destroy(s, fai->hash);
	if (fai->rz) razf_close(fai->rz);
	free(fai);
}

void fai_build(const char *fn)
{
	char *str;
	RAZF *rz;
	FILE *fp;
	faidx_t *fai;
	str = (char*)calloc(strlen(fn) + 5, 1);
	sprintf(str, "%s.fai", fn);
	rz = razf_open(fn, "r");
	assert(rz);
	fai = fai_build_core(rz);
	razf_close(rz);
	fp = fopen(str, "w");
	assert(fp);
	fai_save(fai, fp);
	fclose(fp);
	free(str);
	fai_destroy(fai);
}

faidx_t *fai_load(const char *fn)
{
	char *str;
	FILE *fp;
	faidx_t *fai;
	str = (char*)calloc(strlen(fn) + 5, 1);
	sprintf(str, "%s.fai", fn);
	fp = fopen(str, "r");
	if (fp == 0) {
		fprintf(stderr, "[fai_load] build FASTA index.\n");
		fai_build(fn);
		fp = fopen(str, "r");
		if (fp == 0) {
			free(str);
			return 0;
		}
	}
	fai = fai_read(fp);
	fclose(fp);
	fai->rz = razf_open(fn, "r");
	if (fai->rz == 0) return 0;
	assert(fai->rz);
	free(str);
	return fai;
}

char *fai_fetch(const faidx_t *fai, const char *str, int *len)
{
	char *s, *p, c;
	int i, l, k;
	khiter_t iter;
	faidx1_t val;
	khash_t(s) *h;
	int beg, end;

	beg = end = -1;
	h = fai->hash;
	l = strlen(str);
	p = s = (char*)malloc(l+1);
	/* squeeze out "," */
	for (i = k = 0; i != l; ++i)
		if (str[i] != ',' && !isspace(str[i])) s[k++] = str[i];
	s[k] = 0;
	for (i = 0; i != k; ++i) if (s[i] == ':') break;
	s[i] = 0;
	iter = kh_get(s, h, s); /* get the ref_id */
	if (iter == kh_end(h)) {
		*len = 0;
		free(s); return 0;
	}
	val = kh_value(h, iter);
	if (i == k) { /* dump the whole sequence */
		beg = 0; end = val.len;
	} else {
		for (p = s + i + 1; i != k; ++i) if (s[i] == '-') break;
		beg = atoi(p);
		if (i < k) {
			p = s + i + 1;
			end = atoi(p);
		} else end = val.len;
	}
	if (beg > 0) --beg;
	if (beg >= val.len) beg = val.len;
	if (end >= val.len) end = val.len;
	if (beg > end) beg = end;
	free(s);

	// now retrieve the sequence
	l = 0;
	s = (char*)malloc(end - beg + 2);
	razf_seek(fai->rz, val.offset + beg / val.line_blen * val.line_len + beg % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < end - beg)
		if (isgraph(c)) s[l++] = c;
	s[l] = '\0';
	*len = l;
	return s;
}

int faidx_main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "Usage: faidx <in.fasta> [<reg> [...]]\n");
		return 1;
	} else {
		if (argc == 2) fai_build(argv[1]);
		else {
			int i, j, k, l;
			char *s;
			faidx_t *fai;
			fai = fai_load(argv[1]);
			assert(fai);
			for (i = 2; i != argc; ++i) {
				printf(">%s\n", argv[i]);
				s = fai_fetch(fai, argv[i], &l);
				for (j = 0; j < l; j += 60) {
					for (k = 0; k < 60 && k < l - j; ++k)
						putchar(s[j + k]);
					putchar('\n');
				}
				free(s);
			}
			fai_destroy(fai);
		}
	}
	return 0;
}

#ifdef FAIDX_MAIN
int main(int argc, char *argv[]) { return faidx_main(argc, argv); }
#endif
