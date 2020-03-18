//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_csc_int32: build sparse matrix from int32
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: build sparse matrix from int32 values */

#define SLIP_FREE_ALL                 \
    (*A_handle) = NULL ;                    \
    SLIP_delete_sparse (&A) ;               \
    SLIP_delete_mpz_array(&x_new, nz);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_sparse_csc_int32
(
    SLIP_sparse **A_handle,     // matrix to construct
    int32_t *p,         // The set of column pointers
    int32_t *I,         // set of row indices
    int32_t *x,         // Set of values as doubles
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x and I vectors)
)
{

    SLIP_info info ;
    if (!p || !I || !x || !A_handle)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_sparse *A = slip_create_sparse ( ) ;
    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (A == NULL || x_new == NULL)
    {
        SLIP_FREE_ALL ;
        return (SLIP_OUT_OF_MEMORY) ;
    }

    for (int32_t i = 0; i < nz; i++)
    {
        SLIP_CHECK(SLIP_mpz_set_si(x_new[i], x[i]));
    }
    SLIP_CHECK(SLIP_mpq_set_ui(A->scale, 1, 1));

    SLIP_CHECK(slip_mpz_populate_mat(A, I, p, x_new, n, nz));

    SLIP_delete_mpz_array(&x_new, nz);

    (*A_handle) = A ;
    return SLIP_OK;
}
