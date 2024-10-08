/* $ID: similarity.c, last updated 2020-06-15, F.Osorio */

#include "spatialpack.h"
#include "interface.h"

void
SSIM_coef(double *x, double *y, int *nr, int *nc, double *pars, double *eps, double *stats, double *comp)
{ /* compute structural similarity index for two images */
  int nobs, nrow = *nr, ncol = *nc;
  double xbar, ybar, xvar, yvar, cov;
  double c1, c2, c3;
  double alpha, beta, gamma;
  double luminance, contrast, structure, ssim;

  /* get constants and parameters */
  nobs = nrow * ncol;
  c1   = eps[0];
  c2   = eps[1];
  c3   = eps[2];
  alpha = pars[0];
  beta  = pars[1];
  gamma = pars[2];

  online_covariance(x, y, nobs, &xbar, &ybar, &xvar, &yvar, &cov);

  luminance = (2. * xbar * ybar + c1) / (SQR(xbar) + SQR(ybar) + c1);
  contrast  = (2. * sqrt(xvar) * sqrt(yvar) + c2) / (xvar + yvar + c2);
  structure = (cov + c3) / (sqrt(xvar) * sqrt(yvar) + c3);
  ssim = R_pow(luminance, alpha) * R_pow(contrast, beta) * R_pow(structure, gamma);

  /* save results */
  stats[0] = xbar;
  stats[1] = ybar;
  stats[2] = xvar;
  stats[3] = yvar;
  stats[4] = cov;
  comp[0]  = ssim;
  comp[1]  = luminance;
  comp[2]  = contrast;
  comp[3]  = structure;
}

void
CQ_coef(double *x, double *y, int *nr, int *nc, int *h, double *eps, double *stats, double *comp)
{ /* compute CQ index for two images */
  int nobs, nrow = *nr, ncol = *nc;
  double xbar, ybar, xvar, yvar, cov;
  double c1, c2;
  double luminance, contrast, codispersion;

  /* get constants */
  nobs = nrow * ncol;
  c1   = eps[0];
  c2   = eps[1];

  online_covariance(x, y, nobs, &xbar, &ybar, &xvar, &yvar, &cov);

  luminance = (2. * xbar * ybar + c1) / (SQR(xbar) + SQR(ybar) + c1);
  contrast  = (2. * sqrt(xvar) * sqrt(yvar) + c2) / (xvar + yvar + c2);
  F77_CALL(hcodisp)(x, &nrow, &nrow, &ncol, y, &nrow, h, &codispersion);

  /* save results */
  stats[0] = xbar;
  stats[1] = ybar;
  stats[2] = xvar;
  stats[3] = yvar;
  stats[4] = cov;
  comp[0]  = luminance * contrast * codispersion;
  comp[1]  = luminance;
  comp[2]  = contrast;
  comp[3]  = codispersion;
}
