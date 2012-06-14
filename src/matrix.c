#include "matrix.h"

double
norm_sqr(double *x, int n, int incx)
{   /* sum(x * x) */
    double ans;

    ans = F77_CALL(dnrm2)(&n, x, &incx);
    return R_pow_di(ans, 2);
}

double
dot_product(double *x, int incx, double *y, int incy, int n)
{   /* sum(x * y) */
    double ans;

    ans = F77_CALL(ddot)(&n, x, &incx, y, &incy);
    return ans;
}
