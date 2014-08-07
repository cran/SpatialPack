#ifndef COR_SPATIAL_H
#define COR_SPATIAL_H

#include "base.h"
#include "matrix.h"
#include "stats.h"

/* routines for spatial association */
void cor_tjostheim(double *, double *, double *, double *, DIMS, double, double, double *);
void var_tjostheim(double *, double *, DIMS, double *);
extern void F77_NAME(tjostheim)(double *, double *, double *, double *, int *, double *, double *, double *);

#endif /* COR_SPATIAL_H */
