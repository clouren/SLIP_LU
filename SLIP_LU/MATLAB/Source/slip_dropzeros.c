//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_dropzeros: Drop zeros from a sparse matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"

static bool slip_nonzero (int64_t i, int64_t j, double aij)
{
    return (aij != 0) ;
}

mwIndex slip_dropzeros (mxArray *A)
{
    return (slip_fkeep (A, &slip_nonzero)) ;    /* keep all nonzero entries */
}

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
mwIndex slip_fkeep (mxArray *A, bool (*fkeep) (int64_t, int64_t, double))
{
    if (!A || !mxIsSparse(A)|| !fkeep) return (-1) ;    /* check inputs */

    mwIndex *Ap, *Ai, p, nz = 0;
    mwSize n, j;
    double *Ax ;
    n = mxGetN(A) ;
    Ap = mxGetJc(A) ; Ai = mxGetIr(A) ; Ax = mxGetDoubles(A) ;

    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        Ap [j] = nz ;                       /* record new location of col j */
        for ( ; p < Ap [j+1] ; p++)
        {
            if (fkeep (Ai [p], j, Ax [p]))
            {
                Ax [nz] = Ax [p] ;  /* keep A(i,j) */
                Ai [nz++] = Ai [p] ;
            }
        }
    }
    /* finalize A and remove extra space from A */
    SLIP_realloc(Ai, Ap[n], nz);
    SLIP_realloc(Ax, Ap[n], nz);
    Ap [n] = nz ;
    return (nz) ;
}
