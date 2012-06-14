#ifndef CODISP_COEF_H
#define CODISP_COEF_H

#include "base.h"

/* routines for computation of codispersion coeffcient */
CODISP codisp_init(double *, double *, double *, double *, int *, double *, double *, double *);
void codisp_free(CODISP);
void codisp_coef(double *, double *, DIMS, double *, double *, double *, double *, double *);

#endif /* CODISP_COEF_H */
