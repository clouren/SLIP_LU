#include "SLIP_LU_mex.h"

/* A = A' */
/* This function is modified from CSparse/Source/cs_transpose*/

SLIP_info SLIP_transpose (mxArray *A)
{
    /* check inputs */
    if (!A || !mxIsSparse(A))
    {
        return (SLIP_INCORRECT_INPUT) ;
    }

    mwIndex *Ci, *Ap, *Ai, p, q;
    mwSize m, n, j;
    double *Cx, *Ax ;
    int32_t *Cp, *w;
    m = mxGetM(A) ; n = mxGetN(A) ;
    Ap = mxGetJc(A) ; Ai = mxGetIr(A) ; Ax = mxGetDoubles(A) ;
    w =  (int32_t*) SLIP_calloc (m,    sizeof(int32_t));// get workspace
    Cx = (double *) SLIP_calloc(Ap[n], sizeof(double)) ;// Allocate memory for x
    Cp = (int32_t*) SLIP_calloc(n+1,   sizeof(int32_t));// Initialize p
    Ci = (mwIndex*) SLIP_calloc(Ap[n], sizeof(mwIndex));// Initialize i
    /* out of memory */
    if (!w || !Cx || !Cp || !Ci)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;          /* row counts */
    slip_cumsum (Cp, w, (int32_t) m) ;                     /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            Ci [q = w [Ai [p]]++] = j ; /* place A(i,j) as entry C(j,i) */
            if (Cx) Cx [q] = Ax [p] ;
        }
    }
    SLIP_FREE(Ai);
    mxSetIr(A, Ci);

    SLIP_FREE(Ax);
    mxSetDoubles(A, Cx);

    SLIP_int32_to_mwIndex(Ap, Cp, n+1);
    SLIP_FREE(Cp);

    SLIP_FREE(w);

    return (SLIP_OK) ;
}
