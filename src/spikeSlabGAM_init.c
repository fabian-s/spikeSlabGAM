#define USE_FC_LEN_T
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#ifndef FCONE
# define FCONE
#endif

/* .C calls */
extern void sampler(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, double *, int *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

static const R_CMethodDef CEntries[] = {
    {"sampler", (DL_FUNC) &sampler, 53},
    {NULL, NULL, 0}
};

void R_init_spikeSlabGAM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
