//------------------------------------------------------------------------------
// SLIP_LU/slip_sparse_realloc: double the space for a sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function expands a SLIP LU matrix by doubling its size. This
 * version merely expands x and i and does not initialize/allocate the values!
 */
SLIP_info slip_sparse_realloc
(
    SLIP_sparse* A // the matrix to be expanded
)
{
    if (!A || !A->p || !A->i || !A->x){return SLIP_INCORRECT_INPUT;}
    int32_t nzmax = A->nzmax;
    // Double size of A->x and A->i without initializing the new entries in A->x
    // cannot use SLIP_realloc here, since it frees its input on failure.
    mpz_t *Ax_new = (mpz_t*) SLIP_MEMORY_REALLOC(A->x, 2*nzmax*SIZE_MPZ);

    if (!Ax_new) {return SLIP_OUT_OF_MEMORY;}
    A->x = Ax_new;
    // set newly allocated mpz entries be NULL to avoid potential issue
    for (int32_t i = nzmax; i < 2*nzmax; i++)
    {
        SLIP_MPZ_SET_NULL(A->x[i]);
    }
    A->nzmax = nzmax*2;

    A->i = (int32_t*) SLIP_realloc(A->i, nzmax*sizeof(int32_t),
        2*nzmax*sizeof(int32_t));
    if (!A->i) {return SLIP_OUT_OF_MEMORY;}

    return SLIP_OK;
}
