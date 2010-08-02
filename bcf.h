#ifndef BCF_H
#define BCF_H

#include <stdint.h>
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
	BGZF *fp;
	bcf_hdr_t h;
} bcf_t;

#ifdef __cplusplus
extern "C" {
#endif

	bcf_t *bcf_open(const char *fn, const char *mode);
	int bcf_close(bcf_t *b);
	int bcf_read(bcf_t *bp, bcf1_t *b);
	int bcf_sync(int n_smpl, bcf1_t *b);
	int bcf_write(bcf_t *bp, const bcf1_t *b);
	int bcf_hdr_write(bcf_t *b);
	int bcf_hdr_sync(bcf_hdr_t *b);
	int bcf_destroy(bcf1_t *b);
	char *bcf_fmt(bcf_t *bp, bcf1_t *b);

#ifdef __cplusplus
}
#endif

#endif
