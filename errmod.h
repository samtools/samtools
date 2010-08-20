#ifndef ERRMOD_H
#define ERRMOD_H

#include <stdint.h>

struct __errmod_coef_t;

typedef struct {
	double depcorr;
	struct __errmod_coef_t *coef;
} errmod_t;

errmod_t *errmod_init(float depcorr);
void errmod_destroy(errmod_t *em);
int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q);

#endif
