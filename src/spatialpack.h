#ifndef SPATIALPACK_H
#define SPATIALPACK_H

#include "codisp_coef.h"
#include "cor_spatial.h"
#include "mod_ttest.h"

/* Routines to be called from R */
extern void codisp(double *, double *, double *, double *, int *, double *, double *, double *);
extern void cor_spatial(double *, double *, double *, double *, int *, double *, double *);
extern void modified_ttest(double *, double *, double *, double *, int *, double *, double *, double *, double *, double *);

#endif /* SPATIALPACK_H */
