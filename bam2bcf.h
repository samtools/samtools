#ifndef BAM2BCF_H
#define BAM2BCF_H

#include <stdint.h>
#include "bcf.h"

struct __bcf_callaux_t;
typedef struct __bcf_callaux_t bcf_callaux_t;

typedef struct {
	int depth;
	uint64_t sum_Q2;
	float p[25], esum[5];
} bcf_callret1_t;

typedef struct {
	int a[5]; // alleles: ref, alt, alt2, alt3
	int n, n_alleles, shift, ori_ref, unseen;
	int depth, rmsQ;
	uint8_t *PL;
} bcf_call_t;

#ifdef __cplusplus
extern "C" {
#endif

	bcf_callaux_t *bcf_call_init(double theta);
	void bcf_call_destroy(bcf_callaux_t *bca);
	int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base /*4-bit*/, bcf_callaux_t *bca, bcf_callret1_t *r);
	int bcf_call_combine(int n, const bcf_callret1_t *calls, int ref_base /*4-bit*/, bcf_call_t *call);
	int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b);

#ifdef __cplusplus
}
#endif

#endif
