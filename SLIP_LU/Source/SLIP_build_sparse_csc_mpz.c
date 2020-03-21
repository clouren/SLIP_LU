//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_csc_mpz: build sparse matrix from mpz_t
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose:  build a sparse mpz matrix from another sparse mpz matrix.
 */

#define SLIP_FREE_ALL                 \
    (*A_handle) = NULL ;                    \
    SLIP_delete_sparse (&A) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_sparse_csc_mpz
(
    SLIP_sparse **A_handle,     // matrix to construct
    int32_t *p,         // The set of column pointers
    int32_t *I,         // set of row indices
    mpz_t *x,           // Set of values in full precision integer
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x and I vectors)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_copy (&A, SLIP_CSC, SLIP_MPZ, B, option)
    // instead, and remove this function (or at least make non-user-callable).

    SLIP_info info ;
    if (!p || !I || !x || !A_handle)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------

    SLIP_sparse *A = slip_create_sparse ( ) ;
    if (A == NULL)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    SLIP_CHECK(slip_mpz_populate_mat(A, I, p, x, n, nz));

    SLIP_CHECK(SLIP_mpq_set_ui(A->scale, 1, 1));

    (*A_handle) = A ;
    return SLIP_OK;
}
