#include <R_ext/Rdynload.h>
#include "spatialpack.h"

static const R_CMethodDef CEntries[]  = {
    {"codisp",          (DL_FUNC) &codisp,          8},
    {"cor_spatial",     (DL_FUNC) &cor_spatial,     10},
    {"modified_ttest",  (DL_FUNC) &modified_ttest,  10},
    {NULL, NULL, 0}
};

void R_init_spatools(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
