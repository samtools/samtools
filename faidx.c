#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "faidx.h"
#include "khash.h"

typedef struct {
	int32_t line_len, line_blen;
	int64_t len;
	uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

#ifndef _NO_RAZF
#include "razf.h"
#else
#ifdef _WIN32
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif
#define RAZF FILE
#define razf_read(fp, buf, size) fread(buf, 1, size, fp)
#define razf_open(fn, mode) fopen(fn, mode)
#define razf_close(fp) fclose(fp)
#define razf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define razf_tell(fp) ftello(fp)
#endif
#ifdef _USE_KNETFILE
#include "knetfile.h"
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
	int line_len, line_blen, state;
	int l1, l2;
	faidx_t *idx;
	uint64_t offset;
	int64_t len;

	idx = (faidx_t*)calloc(1, sizeof(faidx_t));
	idx->hash = kh_init(s);
	name = 0; l_name = m_name = 0;
	len = line_len = line_blen = -1; state = 0; l1 = l2 = -1; offset = 0;
	while (razf_read(rz, &c, 1)) {
		if (c == '\n') { // an empty line
			if (state == 1) {
				offset = razf_tell(rz);
				continue;
			} else if ((state == 0 && len < 0) || state == 2) continue;
		}
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
			if (ret == 0) {
				fprintf(stderr, "[fai_build_core] the last entry has no sequence\n");
				free(name); fai_destroy(idx);
				return 0;
			}
			if (c != '\n') while (razf_read(rz, &c, 1) && c != '\n');
			state = 1; len = 0;
			offset = razf_tell(rz);
		} else {
			if (state == 3) {
				fprintf(stderr, "[fai_build_core] inlined empty line is not allowed in sequence '%s'.\n", name);
				free(name); fai_destroy(idx);
				return 0;
			}
			if (state == 2) state = 3;
			l1 = l2 = 0;
			do {
				++l1;
				if (isgraph(c)) ++l2;
			} while ((ret = razf_read(rz, &c, 1)) && c != '\n');
			if (state == 3 && l2) {
				fprintf(stderr, "[fai_build_core] different line length in sequence '%s'.\n", name);
				free(name); fai_destroy(idx);
				return 0;
			}
			++l1; len += l2;
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
#ifdef _WIN32
		fprintf(fp, "%s\t%d\t%ld\t%d\t%d\n", fai->name[i], (int)x.len, (long)x.offset, (int)x.line_blen, (int)x.line_len);
#else
		fprintf(fp, "%s\t%d\t%lld\t%d\t%d\n", fai->name[i], (int)x.len, (long long)x.offset, (int)x.line_blen, (int)x.line_len);
#endif
	}
}

faidx_t *fai_read(FILE *fp)
{
	faidx_t *fai;
	char *buf, *p;
	int len, line_len, line_blen;
#ifdef _WIN32
	long offset;
#else
	long long offset;
#endif
	fai = (faidx_t*)calloc(1, sizeof(faidx_t));
	fai->hash = kh_init(s);
	buf = (char*)calloc(0x10000, 1);
	while (!feof(fp) && fgets(buf, 0x10000, fp)) {
		for (p = buf; *p && isgraph(*p); ++p);
		*p = 0; ++p;
#ifdef _WIN32
		sscanf(p, "%d%ld%d%d", &len, &offset, &line_blen, &line_len);
#else
		sscanf(p, "%d%lld%d%d", &len, &offset, &line_blen, &line_len);
#endif
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

int fai_build(const char *fn)
{
	char *str;
	RAZF *rz;
	FILE *fp;
	faidx_t *fai;
	str = (char*)calloc(strlen(fn) + 5, 1);
	sprintf(str, "%s.fai", fn);
	rz = razf_open(fn, "r");
	if (rz == 0) {
		fprintf(stderr, "[fai_build] fail to open the FASTA file %s\n",fn);
		free(str);
		return -1;
	}
	fai = fai_build_core(rz);
	razf_close(rz);
	fp = fopen(str, "wb");
	if (fp == 0) {
		fprintf(stderr, "[fai_build] fail to write FASTA index %s\n",str);
		fai_destroy(fai); free(str);
		return -1;
	}
	fai_save(fai, fp);
	fclose(fp);
	free(str);
	fai_destroy(fai);
	return 0;
}

#ifdef _USE_KNETFILE
FILE *download_and_open(const char *fn)
{
    const int buf_size = 1 * 1024 * 1024;
    uint8_t *buf;
    FILE *fp;
    knetFile *fp_remote;
    const char *url = fn;
    const char *p;
    int l = strlen(fn);
    for (p = fn + l - 1; p >= fn; --p)
        if (*p == '/') break;
    fn = p + 1;

    // First try to open a local copy
    fp = fopen(fn, "r");
    if (fp)
        return fp;

    // If failed, download from remote and open
    fp_remote = knet_open(url, "rb");
    if (fp_remote == 0) {
        fprintf(stderr, "[download_from_remote] fail to open remote file %s\n",url);
        return NULL;
    }
    if ((fp = fopen(fn, "wb")) == 0) {
        fprintf(stderr, "[download_from_remote] fail to create file in the working directory %s\n",fn);
        knet_close(fp_remote);
        return NULL;
    }
    buf = (uint8_t*)calloc(buf_size, 1);
    while ((l = knet_read(fp_remote, buf, buf_size)) != 0)
        fwrite(buf, 1, l, fp);
    free(buf);
    fclose(fp);
    knet_close(fp_remote);

    return fopen(fn, "r");
}
#endif

faidx_t *fai_load(const char *fn)
{
	char *str;
	FILE *fp;
	faidx_t *fai;
	str = (char*)calloc(strlen(fn) + 5, 1);
	sprintf(str, "%s.fai", fn);

#ifdef _USE_KNETFILE
    if (strstr(fn, "ftp://") == fn || strstr(fn, "http://") == fn)
    {
        fp = download_and_open(str);
        if ( !fp )
        {
            fprintf(stderr, "[fai_load] failed to open remote FASTA index %s\n", str);
            free(str);
            return 0;
        }
    }
    else
#endif
        fp = fopen(str, "rb");
	if (fp == 0) {
		fprintf(stderr, "[fai_load] build FASTA index.\n");
		fai_build(fn);
		fp = fopen(str, "rb");
		if (fp == 0) {
			fprintf(stderr, "[fai_load] fail to open FASTA index.\n");
			free(str);
			return 0;
		}
	}

	fai = fai_read(fp);
	fclose(fp);

	fai->rz = razf_open(fn, "rb");
	free(str);
	if (fai->rz == 0) {
		fprintf(stderr, "[fai_load] fail to open FASTA file.\n");
		return 0;
	}
	return fai;
}

char *fai_fetch(const faidx_t *fai, const char *str, int *len)
{
	char *s, c;
	int i, l, k, name_end;
	khiter_t iter;
	faidx1_t val;
	khash_t(s) *h;
	int beg, end;

	beg = end = -1;
	h = fai->hash;
	name_end = l = strlen(str);
	s = (char*)malloc(l+1);
	// remove space
	for (i = k = 0; i < l; ++i)
		if (!isspace(str[i])) s[k++] = str[i];
	s[k] = 0; l = k;
	// determine the sequence name
	for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
	if (i >= 0) name_end = i;
	if (name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for (i = name_end + 1; i < l; ++i) {
			if (s[i] == '-') ++n_hyphen;
			else if (!isdigit(s[i]) && s[i] != ',') break;
		}
		if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
		s[name_end] = 0;
		iter = kh_get(s, h, s);
		if (iter == kh_end(h)) { // cannot find the sequence name
			iter = kh_get(s, h, str); // try str as the name
			if (iter == kh_end(h)) {
				*len = 0;
			free(s); return 0;
			} else s[name_end] = ':', name_end = l;
		}
	} else iter = kh_get(s, h, str);
	val = kh_value(h, iter);
	// parse the interval
	if (name_end < l) {
		for (i = k = name_end + 1; i < l; ++i)
			if (s[i] != ',') s[k++] = s[i];
		s[k] = 0;
		beg = atoi(s + name_end + 1);
		for (i = name_end + 1; i != k; ++i) if (s[i] == '-') break;
		end = i < k? atoi(s + i + 1) : val.len;
		if (beg > 0) --beg;
	} else beg = 0, end = val.len;
	if (beg >= val.len) beg = val.len;
	if (end >= val.len) end = val.len;
	if (beg > end) beg = end;
	free(s);

	// now retrieve the sequence
	l = 0;
	s = (char*)malloc(end - beg + 2);
	razf_seek(fai->rz, val.offset + beg / val.line_blen * val.line_len + beg % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < end - beg && !fai->rz->z_err)
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
			if (fai == 0) return 1;
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

int faidx_fetch_nseq(const faidx_t *fai) 
{
	return fai->n;
}

char *faidx_fetch_seq(const faidx_t *fai, char *c_name, int p_beg_i, int p_end_i, int *len)
{
	int l;
	char c;
    khiter_t iter;
    faidx1_t val;
	char *seq=NULL;

    // Adjust position
    iter = kh_get(s, fai->hash, c_name);
    if(iter == kh_end(fai->hash)) return 0;
    val = kh_value(fai->hash, iter);
	if(p_end_i < p_beg_i) p_beg_i = p_end_i;
    if(p_beg_i < 0) p_beg_i = 0;
    else if(val.len <= p_beg_i) p_beg_i = val.len - 1;
    if(p_end_i < 0) p_end_i = 0;
    else if(val.len <= p_end_i) p_end_i = val.len - 1;

    // Now retrieve the sequence 
	l = 0;
	seq = (char*)malloc(p_end_i - p_beg_i + 2);
	razf_seek(fai->rz, val.offset + p_beg_i / val.line_blen * val.line_len + p_beg_i % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < p_end_i - p_beg_i + 1)
		if (isgraph(c)) seq[l++] = c;
	seq[l] = '\0';
	*len = l;
	return seq;
}

#ifdef FAIDX_MAIN
int main(int argc, char *argv[]) { return faidx_main(argc, argv); }
#endif
