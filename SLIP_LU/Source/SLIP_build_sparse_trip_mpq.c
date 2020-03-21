//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_trip_mpq: build sparse matrix from mpq_t
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose:  build sparse mpz matrix from mpq triplets
 */

#define SLIP_FREE_ALL                 \
    (*A_handle) = NULL ;                    \
    SLIP_delete_sparse (&A) ;               \
    SLIP_delete_mpz_array(&x_new, nz);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_sparse_trip_mpq
(
    SLIP_sparse **A_handle,     // matrix to construct
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    mpq_t *x,           // Set of values as rational numbers
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x, I, and J vectors)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_copy (&C, SLIP_CSC, SLIP_MPZ, B, option)
    // instead, and remove this function (or at least make non-user-callable).

    SLIP_info info ;
    if (!I || !J || !A_handle || !x || n <= 0 || nz <= 0)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------

    SLIP_sparse *A = slip_create_sparse ( ) ;
    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (A == NULL || x_new == NULL)
    {
        SLIP_FREE_ALL ;
        return (SLIP_OUT_OF_MEMORY) ;
    }

    SLIP_CHECK(slip_expand_mpq_array(x_new, x, A->scale, nz));

    SLIP_CHECK(slip_trip_to_mat(A, I, J, x_new, n, nz));

    SLIP_delete_mpz_array(&x_new, nz);

    (*A_handle) = A ;
    return SLIP_OK;
}

