#include "codisp_coef.h"

void
codisp(double *x, double *y, double *xpos, double *ypos, int *pdims, 
    double *upper_bounds, double *card, double *coef)
{
    CODISP obj;

    obj = codisp_init(x, y, xpos, ypos, pdims, upper_bounds, card, coef);
    codisp_coef(obj->data->x, obj->data->y, obj->data->dims, obj->data->xpos,
                obj->data->ypos, obj->data->upper_bounds, obj->data->card,
                obj->coef);
    codisp_free(obj);
}

CODISP
codisp_init(double *x, double *y, double *xpos, double *ypos, int *pdims,
    double *upper_bounds, double *card, double *coef)
{   /* data object */
    CODISP obj;
    int do_half = 1;

    obj = (CODISP) Calloc(1, CODISP_struct);
    obj->data = data_init(x, y, xpos, ypos, pdims, do_half, upper_bounds, card);
    obj->coef = coef;
    return obj;
}

void
codisp_free(CODISP this)
{   /* destructor for a codisp object */
    data_free(this->data);
    Free(this);
}

void
codisp_coef(double *x, double *y, DIMS dims, double *xpos, double *ypos,
    double *upper_bounds, double *card, double *coef)
{
    int i, j, k, which_class;
    double accum, distance, dx, dy, sxx, syy, sxy;
    
    for (k = 0; k < dims->nclass; k++) {
        accum = 0.0;
        sxx = syy = sxy = 0.0;
        for (i = 0; i < dims->n; i++) {
            for (j = i + 1; j < dims->n; j++) {
                dx = (xpos[i] - xpos[j]);
                dy = (ypos[i] - ypos[j]);
                distance = hypot(dx, dy);
                which_class = find_interval(upper_bounds, dims->nclass, distance);
                if (which_class == k) {
                    accum += 1.0;
                    dx = (x[i] - x[j]);
                    dy = (y[i] - y[j]);
                    sxx += SQR(dx);
                    syy += SQR(dy);
                    sxy += dx * dy;
                }
            }
        }
        sxx /= accum;
        syy /= accum;
        sxy /= accum;
        card[k] = accum;
        if (accum)
            coef[k] = sxy / sqrt(sxx * syy);
    }
}
