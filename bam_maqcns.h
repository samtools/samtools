#ifndef BAM_MAQCNS_H
#define BAM_MAQCNS_H

#include "glf.h"

struct __bmc_aux_t;

typedef struct {
	float het_rate, theta;
	int n_hap, cap_mapQ, is_soap;

	float eta, q_r;
	double *fk, *coef;
	double *lhet;
	struct __bmc_aux_t *aux;
} bam_maqcns_t;

typedef struct {
	int q_indel;
	float r_indel;
	// hidden parameters, unchangeable from command line
	int mm_penalty, indel_err, ambi_thres;
} bam_maqindel_opt_t;

typedef struct {
	int indel1, indel2;
	int cnt1, cnt2, cnt_anti;
	int cnt_ref, cnt_ambi;
	char *s[2];
	//
	int gt, gl[2];
	int q_cns, q_ref;
} bam_maqindel_ret_t;

#ifdef __cplusplus
extern "C" {
#endif

	bam_maqcns_t *bam_maqcns_init();
	void bam_maqcns_prepare(bam_maqcns_t *bm);
	void bam_maqcns_destroy(bam_maqcns_t *bm);
	glf1_t *bam_maqcns_glfgen(int n, const bam_pileup1_t *pl, uint8_t ref_base, bam_maqcns_t *bm);
	uint32_t bam_maqcns_call(int n, const bam_pileup1_t *pl, bam_maqcns_t *bm);
	// return: cns<<28 | cns2<<24 | mapQ<<16 | cnsQ<<8 | cnsQ2
	uint32_t glf2cns(const glf1_t *g, int q_r);

	bam_maqindel_opt_t *bam_maqindel_opt_init();
	bam_maqindel_ret_t *bam_maqindel(int n, int pos, const bam_maqindel_opt_t *mi, const bam_pileup1_t *pl, const char *ref,
									 int _n_types, int *_types);
	void bam_maqindel_ret_destroy(bam_maqindel_ret_t*);

#ifdef __cplusplus
}
#endif

#endif
