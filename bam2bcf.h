#ifndef BAM2BCF_H
#define BAM2BCF_H

#include <stdint.h>
#include "errmod.h"
#include "bcftools/bcf.h"

#define BAM2BCF_INDELREG_THRES 13

typedef struct __bcf_callaux_t {
	int max_bases, capQ, min_baseQ;
	int openQ, extQ, tandemQ;
	// for internal uses
	uint16_t *bases;
	errmod_t *e;
} bcf_callaux_t;

typedef struct {
	int depth, indelreg;
	float ins_len, del_len;
	int qsum[4];
	int anno[16];
	float p[25];
} bcf_callret1_t;

typedef struct {
	int a[5]; // alleles: ref, alt, alt2, alt3
	int n, n_alleles, shift, ori_ref, unseen;
	int indelreg;
	float ins_len, del_len;
	int depth, anno[16];
	uint8_t *PL;
} bcf_call_t;

#ifdef __cplusplus
extern "C" {
#endif

	bcf_callaux_t *bcf_call_init(double theta, int min_baseQ);
	void bcf_call_destroy(bcf_callaux_t *bca);
	int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base /*4-bit*/, bcf_callaux_t *bca, bcf_callret1_t *r);
	int bcf_call_combine(int n, const bcf_callret1_t *calls, int ref_base /*4-bit*/, bcf_call_t *call);
	int bcf_call_combine_gap(int n, const bcf_callret1_t *calls, bcf_call_t *call);
	int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b, bcf_callret1_t *bcr, int is_SP);
	int bcf_call_glfgen_gap(int pos, int _n, const bam_pileup1_t *pl, bcf_callaux_t *bca, bcf_callret1_t *r);

#ifdef __cplusplus
}
#endif

#endif
