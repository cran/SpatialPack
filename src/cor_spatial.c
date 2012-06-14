#include "cor_spatial.h"

void
cor_spatial(double *x, double *y, double *xpos, double *ypos, int *pdims,
    double *cor, double *var)
{
    DIMS dims;

    dims = dimension(pdims);
    cor_tjostheim(x, y, dims, cor);
    var_tjostheim(xpos, ypos, dims, var);
    dimension_free(dims);
}

void
cor_tjostheim(double *x, double *y, DIMS dims, double *cor)
{
    int i;
    double ans = 0.0, sx = 0.0, sy = 0.0;

    for (i = 0; i < dims->n; i++) {
        sx += SQR(x[i]) + SQR(y[i]);
        sy += SQR(x[i + dims->n]) + SQR(y[i + dims->n]);
    }
    ans += dot_product(x, 1, x + dims->n, 1, dims->n);
    ans += dot_product(y, 1, y + dims->n, 1, dims->n);
    ans /= sqrt(sx * sy);
    *cor = ans;
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
