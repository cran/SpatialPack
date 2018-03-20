#include "similarity.h"

void
SSIM(double *x, double *y, double *params, double *eps, double *stats, double *comp)
{
  int nobs;
  double xbar, ybar, xvar, yvar, cov;
  double c1, c2, c3;
  double alpha, beta, gamma;
  double luminance, contrast, corr, ssim;

  /* get constants and parameters */
  nobs = (int) eps[0];
  c1   = eps[1];
  c2   = eps[2];
  c3   = eps[3];
  alpha = params[0];
  beta  = params[1];
  gamma = params[2];

  online_covariance(x, y, nobs, &xbar, &ybar, &xvar, &yvar, &cov);

  luminance = (2. * xbar * ybar + c1) / (SQR(xbar) + SQR(ybar) + c1);
  contrast  = (2. * sqrt(xvar) * sqrt(yvar) + c2) / (xvar + yvar + c2);
  corr = (cov + c3) / (sqrt(xvar) * sqrt(yvar) + c3);
  ssim = luminance * contrast * corr;

  /* save results */
  stats[0] = xbar;
  stats[1] = ybar;
  stats[2] = xvar;
  stats[3] = yvar;
  stats[4] = cov;
  comp[0]  = ssim;
  comp[1]  = luminance;
  comp[2]  = contrast;
  comp[3]  = corr;
}
