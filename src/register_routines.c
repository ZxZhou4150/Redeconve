#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
// #include <stdlib.h> // for NULL

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

// Fortran functions
extern void F77_NAME(Hessian_f90)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"Hessian_f", (DL_FUNC)&F77_NAME(Hessian_f90), 5},
    {NULL, NULL, 0}};

void R_init_Redeconve(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
