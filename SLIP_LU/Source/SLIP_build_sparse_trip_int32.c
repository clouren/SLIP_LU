//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_trip_int32: build sparse matrix from int32
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose:  build a sparse mpz matrix from int32 triplets
 */

#define SLIP_FREE_ALL                       \
    (*A_handle) = NULL ;                    \
    SLIP_delete_sparse (&A) ;               \
    SLIP_delete_mpz_array(&x_new, nz);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_sparse_trip_int32
(
    SLIP_sparse **A_handle,     // matrix to construct
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    int32_t *x,         // Set of values in int32
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

    for (int32_t k = 0; k < nz; k++)
    {
        SLIP_CHECK(SLIP_mpz_set_si(x_new[k], x[k]));
    }
    SLIP_CHECK(SLIP_mpq_set_ui(A->scale, 1,1));

    SLIP_CHECK(slip_trip_to_mat(A, I, J, x_new, n, nz));

    SLIP_delete_mpz_array(&x_new, nz);

    (*A_handle) = A ;
    return SLIP_OK;
}
