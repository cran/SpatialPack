#ifndef COR_SPATIAL_H
#define COR_SPATIAL_H

#include "base.h"
#include "matrix.h"
#include "stats.h"

/* routines for spatial association */
void cor_tjostheim(double *, double *, DIMS, double *);
void var_tjostheim(double *, double *, DIMS, double *);

#endif /* COR_SPATIAL_H */
