//------------------------------------------------------------------------------
// SLIP_LU/slip_sparse_collapse: shrink space required by a sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function collapses a SLIP matrix. Essentially it shrinks the
 * size of x and i. so that they only take up the number of elements in the
 * matrix. For example if A->nzmax = 1000 but A->nz = 500, A->i and A->x are of size
 * 1000, so this function shrinks them to size 500
 */
SLIP_info slip_sparse_collapse
(
    SLIP_sparse* A // matrix to be shrunk
)
{
    if (!A || !A->p || !A->i || !A->x) {return SLIP_INCORRECT_INPUT;}
    int32_t nz = A->nz;
    // Shrink A->i and A->x such that they're of size nz.
    // SLIP_realloc cannot fail in this case since the space is shrinking.
    A->i = (int32_t*) SLIP_realloc(A->i, A->nzmax*sizeof(int32_t),
        nz*sizeof(int32_t));
    // SLIP_realloc can be used here, since it would not fail when shrinking.
    A->x = (mpz_t*) SLIP_realloc(A->x, A->nzmax*SIZE_MPZ, nz*SIZE_MPZ);
    A->nzmax = nz;
    return SLIP_OK;
}
