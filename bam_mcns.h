#ifndef BAM_MCNS_H
#define BAM_MCNS_H

#include "bam.h"

struct __mc_aux_t;
typedef struct __mc_aux_t mc_aux_t;

#ifdef __cplusplus
extern "C" {
#endif

	mc_aux_t *mc_init(int n);
	void mc_destroy(mc_aux_t *ma);
	double mc_freq0(int ref, int *n, const bam_pileup1_t **plp, mc_aux_t *ma, int *_ref, int *alt);
	double mc_freq_iter(double f0, mc_aux_t *ma);
	double mc_ref_prob(mc_aux_t *ma);
	int mc_call_gt(const mc_aux_t *ma, double f0, int k);

#ifdef __cplusplus
}
#endif

#endif
