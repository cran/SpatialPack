/* $ID: stats.h, last updated 2020-06-08, F.Osorio */

#ifndef STATS_H
#define STATS_H

#include "base.h"

/* routines for the computations of sample statistics */
extern void mean_and_var(double *, int, double *, double *);
extern void online_covariance(double *, double *, int, double *, double *, double *, double *, double *);
extern void F77_SUB(moments)(double *, int *, double *, double *);
extern double F77_SUB(median)(double *, int *);

#endif /* STATS_H */
