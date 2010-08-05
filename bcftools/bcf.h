#ifndef BCF_H
#define BCF_H

#include <stdint.h>
#include <zlib.h>
#include "bgzf.h"

typedef struct {
	int fmt, len; // len is the unit length
	void *data;
	// derived info: fmt, len
} bcf_ginfo_t;

typedef struct {
	int32_t tid, pos;
	uint32_t qual:8, l_str:24;
	int m_str;
	char *str, *ref, *alt, *flt, *info, *fmt; // fmt, ref, alt and info point to str
	int n_gi, m_gi;
	bcf_ginfo_t *gi;
	int n_alleles;
	// derived info: ref, alt, flt, info, fmt, n_gi, n_alleles
} bcf1_t;

typedef struct {
	int32_t n_ref, n_smpl;
	int32_t l_nm;
	int32_t l_smpl;
	int32_t l_txt;
	char *name, *sname, *txt;
	char **ns, **sns;
	// derived info: n_ref, n_smpl, ns, sns
} bcf_hdr_t;

typedef struct {
	int is_vcf;
	void *v;
	BGZF *fp;
} bcf_t;

struct __bcf_idx_t;
typedef struct __bcf_idx_t bcf_idx_t;

#ifdef __cplusplus
extern "C" {
#endif

	bcf_t *bcf_open(const char *fn, const char *mode);
	int bcf_close(bcf_t *b);
	int bcf_read(bcf_t *bp, const bcf_hdr_t *h, bcf1_t *b);
	int bcf_sync(int n_smpl, bcf1_t *b);
	int bcf_write(bcf_t *bp, const bcf_hdr_t *h, const bcf1_t *b);
	bcf_hdr_t *bcf_hdr_read(bcf_t *b);
	int bcf_hdr_write(bcf_t *b, const bcf_hdr_t *h);
	int bcf_hdr_sync(bcf_hdr_t *b);
	void bcf_hdr_destroy(bcf_hdr_t *h);
	int bcf_destroy(bcf1_t *b);
	char *bcf_fmt(const bcf_hdr_t *h, bcf1_t *b);

	int vcf_close(bcf_t *bp);

	void *bcf_build_refhash(bcf_hdr_t *h);
	void bcf_str2id_destroy(void *_hash);
	int bcf_str2id(void *_hash, const char *str);

	int bcf_idx_build(const char *fn);
	uint64_t bcf_idx_query(const bcf_idx_t *idx, int tid, int beg, int end);
	int bcf_parse_region(void *str2id, const char *str, int *tid, int *begin, int *end);
	bcf_idx_t *bcf_idx_load(const char *fn);
	void bcf_idx_destroy(bcf_idx_t *idx);

#ifdef __cplusplus
}
#endif

#endif
