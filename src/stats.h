#ifndef SPATOOLS_STATS_H
#define SPATOOLS_STATS_H

#include "base.h"

/* routines for the computations of sample statistics */
void mean_and_var(double *, int, double *, double *);
void online_covariance(double *, double *, int, double *, double *, double *, double *, double *);

#endif /* SPATOOLS_STATS_H */
