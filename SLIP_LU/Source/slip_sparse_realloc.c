//------------------------------------------------------------------------------
// SLIP_LU/slip_sparse_realloc: double the space for a sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function expands a CSC SLIP_matrix by doubling its size. This
 * version merely expands x and i and does not initialize/allocate the values!
 * The only purpose of this function is for the factorization, it does not work 
 * for general sparse matrices
 */

#include "SLIP_LU_internal.h"

SLIP_info slip_sparse_realloc
(
    SLIP_matrix* A // the matrix to be expanded
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    
    int64_t nzmax = A->nzmax;
    mpz_t *Ax_new = NULL;

    // Double size of A->x and A->i without initializing the new entries in A->x
    // cannot use SLIP_realloc here, since it frees its input on failure.
    Ax_new = (mpz_t*) SLIP_MEMORY_REALLOC (A->x.mpz, 2*nzmax*sizeof(mpz_t));
    if (!Ax_new)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    A->x.mpz = Ax_new;

    // set newly allocated mpz entries to NULL
    for (int64_t i = nzmax; i < 2*nzmax; i++)
    {
        SLIP_MPZ_SET_NULL (A->x.mpz[i]) ;
    }

    A->i = (int64_t*) SLIP_realloc (A->i, nzmax*sizeof(int64_t),
        2*nzmax*sizeof(int64_t)) ;
    if (!A->i)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    A->nzmax = nzmax*2;
    return (SLIP_OK) ;
}

