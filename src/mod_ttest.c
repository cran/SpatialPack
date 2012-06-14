#include "mod_ttest.h"

void
modified_ttest(double *x, double *y, double *xpos, double *ypos, int *pdims,
    double *cor, double *upper_bounds, double *card, double *imoran, double *stats)
{
    TTEST obj;

    obj = mod_ttest_init(x, y, xpos, ypos, pdims, cor, upper_bounds, card,
                         imoran, stats);
    mod_ttest(obj->data->x, obj->data->y, obj->data->dims, obj->data->xpos,
              obj->data->ypos, obj->data->upper_bounds, obj->cor, obj->data->card,
              obj->imoran, obj->stats);
    mod_ttest_free(obj);
}

TTEST
mod_ttest_init(double *x, double *y, double *xpos, double *ypos, int *pdims,
    double *cor, double *upper_bounds, double *card, double *imoran, double *stats)
{   /* modified t-test object */
    TTEST obj;
    int do_half = 0;

    obj = (TTEST) Calloc(1, TTEST_struct);
    obj->data = data_init(x, y, xpos, ypos, pdims, do_half, upper_bounds, card);
    obj->cor = cor;
    obj->imoran = imoran;
    obj->stats = stats;
    return obj;
}

void
mod_ttest_free(TTEST this)
{   /* destructor for a t-test object */
    data_free(this->data);
    Free(this);
}

void
MoranI(double *x, double *y, DIMS dims, double *xpos, double *ypos,
    double *upper_bounds, double *card, double *index)
{   /* Moran's I */
    int i, j, k, which_class;
    double dx, dy, distance, sx, sy, xbar, xvar, ybar, yvar, wts;

    mean_and_var(x, dims->n, &xbar, &xvar);
    mean_and_var(y, dims->n, &ybar, &yvar);

    for (k = 0; k < dims->nclass; k++) {
        sx = 0.0;
        sy = 0.0;
        wts = 0.0;
        for (j = 0; j < dims->n; j++) {
            for (i = j + 1; i < dims->n; i++) {
                dx = (xpos[i] - xpos[j]);
                dy = (ypos[i] - ypos[j]);
                distance = hypot(dx, dy);
                which_class = find_interval(upper_bounds, dims->nclass, distance);
                if (which_class == k) {
                    wts += 1.0;
                    sx += (x[i] - xbar) * (x[j] - xbar);
                    sy += (y[i] - ybar) * (y[j] - ybar);
                }
            }
        }
        index[k] = (sx / wts) / xvar;
        index[k + dims->nclass] = (sy / wts) / yvar;
        card[k] = wts;
    }
}

double
corrected_df(double *xpos, double *ypos, DIMS dims, double *upper_bounds,
    double *imoran)
{
    int i, j, pos;
    double corx, cory, dx, dy, distance, rxx, ryy, sxx, syy, sxy, trxx, tryy, trxy, ans;
    
    /* initialization of correlation matrices */
    sxx = syy = sxy = trxy = 0.0;
    for (j = 0; j < dims->n; j++) {
        rxx = ryy = 0.0;
        for (i = 0; i < dims->n; i++) {
            if (i != j) {
                dx = (xpos[i] - xpos[j]);
                dy = (ypos[i] - ypos[j]);
                distance = hypot(dx, dy);
                pos = find_interval(upper_bounds, dims->nclass, distance);
                corx = imoran[pos];
                cory = imoran[pos + dims->nclass];
            }
            else
                corx = cory = 1.0;
            rxx  += corx;
            ryy  += cory;
            trxy += corx * cory;
        }
        sxx += rxx;
        syy += ryy;
        sxy += rxx * ryy;
    }
    
    /* computation of the corrected number of degrees of freedom */
    trxx = (double) dims->n - sxx / dims->n;
    tryy = (double) dims->n - syy / dims->n;
    trxy += (sxx * syy / dims->n - 2.0 * sxy) / dims->n;
    ans = trxx * tryy / trxy - 1.0;

    return ans;
}

void
mod_ttest(double *x, double *y, DIMS dims, double *xpos, double *ypos, 
    double *upper_bounds, double *cor, double *card, double *imoran, double *stats)
{
    int lower_tail = 0, log_p = 0;
    double R, F, df, pval;

    MoranI(x, y, dims, xpos, ypos, upper_bounds, card, imoran);

    df = corrected_df(xpos, ypos, dims, upper_bounds, imoran);
    R  = *cor; /* only for p = 2! */
    F  = df * SQR(R) / (1.0 - SQR(R));
    pval = pf(F, 1.0, df, lower_tail, log_p);

    /* save results */
    stats[0] = F;
    stats[1] = df;
    stats[2] = pval;
}
