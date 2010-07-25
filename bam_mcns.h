#ifndef BAM_MCNS_H
#define BAM_MCNS_H

#include "bam.h"

struct __mc_aux_t;
typedef struct __mc_aux_t mc_aux_t;

typedef struct {
	// O(n)
	int ref, alt, alt2;
	double f_em, f_naive, f_nielsen;
	// O(n^2)
	double PD, p_ref, f_exp;
	// O(n^3)
	double f_map, p_map; // map=maximum a posterior
} mc_rst_t;

#define MC_PTYPE_FULL  1
#define MC_PTYPE_COND2 2
#define MC_PTYPE_FLAT  3

#ifdef __cplusplus
extern "C" {
#endif

	mc_aux_t *mc_init(int n);
	void mc_init_prior(mc_aux_t *ma, int type, double theta);
	void mc_destroy(mc_aux_t *ma);
	int mc_cal(int ref, int *n, const bam_pileup1_t **plp, mc_aux_t *ma, mc_rst_t *rst, int level);
	int mc_call_gt(const mc_aux_t *ma, double f0, int k);
	void mc_dump_afs(mc_aux_t *ma);

#ifdef __cplusplus
}
#endif

#endif
