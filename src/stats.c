#include "stats.h"

void
mean_and_var(double *x, int nobs, double *mean, double *var)
{   /* computes the sample mean and variance using an online algorithm */
    int i, n = 0;
    double accum = 0.0, diff;
    
    *mean = 0.0;
    for (i = 0; i < nobs; i++) {
        n++;
        diff = x[i] - *mean;
        *mean += diff / n;
        accum += diff * (x[i] - *mean);
    }
    *var = accum / n;
}
