#ifndef MOD_TTEST_H
#define MOD_TTEST_H

#include "base.h"
#include "matrix.h"
#include "stats.h"

/* routines for computation of modified t-test */
TTEST mod_ttest_init(double *, double *, double *, double *, int *, double *, double *, double *, double *, double *);
void mod_ttest_free(TTEST);
void MoranI(double *, double *, DIMS, double *, double *, double *, double *, double *);
double corrected_df(double *, double *, DIMS, double *, double *);
void mod_ttest(double *, double *, DIMS, double *, double *, double *, double *, double *, double *, double *);

#endif /* MOD_TTEST_H */
