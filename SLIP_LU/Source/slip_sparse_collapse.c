//------------------------------------------------------------------------------
// SLIP_LU/slip_sparse_collapse: shrink space required by a CSC mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function collapses a SLIP matrix. Essentially it shrinks the
 * size of x and i. so that they only take up the number of elements in the
 * matrix. For example if A->nzmax = 1000 but A->nz = 500, A->i and A->x are of
 * size 1000, so this function shrinks them to size 500.
 * This is only valid in the factorization routines for sparse csc mpz matrices
 */

#include "SLIP_LU_internal.h"
SLIP_info slip_sparse_collapse
(
    SLIP_matrix* A // matrix to be shrunk
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);

    //--------------------------------------------------------------------------

    int64_t nz = A->p[n];
    // Shrink A->i and A->x such that they're of size nz.  These calls to
    // SLIP_realloc cannot fail since the space is shrinking.
    A->i = (int64_t*) SLIP_realloc(A->i, A->nzmax*sizeof(int64_t),
        nz*sizeof(int64_t));
    A->x.mpz = (mpz_t*) SLIP_realloc(A->x.mpz, A->nzmax*sizeof(mpz_t),
        nz*sizeof(mpz_t));
    A->nzmax = nz;
    return SLIP_OK;
}
