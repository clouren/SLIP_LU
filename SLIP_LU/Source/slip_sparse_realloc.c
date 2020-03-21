//------------------------------------------------------------------------------
// SLIP_LU/slip_sparse_realloc: double the space for a sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function expands a SLIP LU matrix by doubling its size. This
 * version merely expands x and i and does not initialize/allocate the values!
 */

// TODO: rename this SLIP_matrix_reallocate, and extend to all matrix kinds
// and types?

#include "SLIP_LU_internal.h"

SLIP_info slip_sparse_realloc
(
    SLIP_sparse* A // the matrix to be expanded
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO could extend this to any data type, but it's only needed for
    // L and U, in SLIP_LU_factorize.

    ASSERT (A != NULL && A->kind = SLIP_CSC && A->type == SLIP_MPZ) ;

    if (!A || !A->p || !A->i || !A->x){return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------

    int64_t nzmax = A->nzmax;

    // Double size of A->x and A->i without initializing the new entries in A->x
    // cannot use SLIP_realloc here, since it frees its input on failure.
    mpz_t *Ax_new = (mpz_t*) SLIP_MEMORY_REALLOC (A->x, 2*nzmax*sizeof(mpz_t));
    if (!Ax_new)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    A->x = Ax_new;

    // set newly allocated mpz entries to NULL
    for (int64_t i = nzmax; i < 2*nzmax; i++)
    {
        SLIP_MPZ_SET_NULL (A->x[i]) ;
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

