#include "cor_spatial.h"

void
cor_spatial(double *xposx, double *xposy, double *yposx, double *yposy, double *bars,
            double *xpos, double *ypos, int *pdims, double *cor, double *var)
{
    DIMS dims;
    double xbar, ybar;

    dims = dimension(pdims);
    xbar = bars[0];
    ybar = bars[1];
    cor_tjostheim(xposx, xposy, yposx, yposy, dims, xbar, ybar, cor);
    var_tjostheim(xpos, ypos, dims, var);
    dimension_free(dims);
}

void
cor_tjostheim(double *xposx, double *xposy, double *yposx, double *yposy, DIMS dims,
              double xbar, double ybar, double *cor)
{   
    F77_CALL(tjostheim)(xposx, xposy, yposx, yposy, &(dims->n), &xbar, &ybar, cor);
}

void
var_tjostheim(double *xpos, double *ypos, DIMS dims, double *var)
{
    double ans = 0.0, sxx, syy, sxy;

    sxx = norm_sqr(xpos, dims->n, 1);
    syy = norm_sqr(ypos, dims->n, 1);
    sxy = dot_product(xpos, 1, ypos, 1, dims->n);
    ans += SQR(sxx) + SQR(syy) + 2.0 * SQR(sxy);
    ans /= (dims->n - 1.0) * SQR(sxx + syy);
    *var = ans;
}
