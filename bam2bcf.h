#ifndef BAM2BCF_H
#define BAM2BCF_H

#include <stdint.h>
#include "errmod.h"
#include "bcftools/bcf.h"

#define B2B_INDEL_NULL 10000

typedef struct __bcf_callaux_t {
	int capQ, min_baseQ;
	int openQ, extQ, tandemQ; // for indels
	int min_support; // for collecting indel candidates
	double min_frac; // for collecting indel candidates
	// for internal uses
	int max_bases;
	int indel_types[4];
	int maxins, indelreg;
	char *inscns;
	uint16_t *bases;
	errmod_t *e;
	void *rghash;
} bcf_callaux_t;

typedef struct {
	int depth, ori_depth, qsum[4];
	int anno[16];
	float p[25];
} bcf_callret1_t;

typedef struct {
	int a[5]; // alleles: ref, alt, alt2, alt3
	int n, n_alleles, shift, ori_ref, unseen;
	int anno[16], depth, ori_depth;
	uint8_t *PL;
} bcf_call_t;

#ifdef __cplusplus
extern "C" {
#endif

	bcf_callaux_t *bcf_call_init(double theta, int min_baseQ);
	void bcf_call_destroy(bcf_callaux_t *bca);
	int bcf_call_glfgen(int _n, const bam_pileup1_t *pl, int ref_base, bcf_callaux_t *bca, bcf_callret1_t *r);
	int bcf_call_combine(int n, const bcf_callret1_t *calls, int ref_base /*4-bit*/, bcf_call_t *call);
	int bcf_call2bcf(int tid, int pos, bcf_call_t *bc, bcf1_t *b, bcf_callret1_t *bcr, int is_SP,
					 const bcf_callaux_t *bca, const char *ref);
	int bcf_call_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref,
						  const void *rghash);

#ifdef __cplusplus
}
#endif

#endif
