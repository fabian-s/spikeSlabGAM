#define USE_FC_LEN_T
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#ifndef FCONE
# define FCONE
#endif

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void sampler(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"sampler", (DL_FUNC) &sampler, 53},
    {NULL, NULL, 0}
};

void R_init_spikeSlabGAM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
