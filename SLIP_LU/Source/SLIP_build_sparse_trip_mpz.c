//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_trip_mpz: build sparse matrix from mpz_t
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function will allow the user to take a matrix of their defined
 * type (in this case mpz_t) and convert it from their triplet form to our data
 * structure. The integrity of the user defined arrays are maintained
 * (therefore, one would need to delete these arrays).
 *
 * On output, the SLIP_sparse* A contains the user's matrix.
 */

#define SLIP_FREE_WORKSPACE                 \
    (*A_handle) = NULL ;                    \
    SLIP_delete_sparse (&A) ;               \

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_sparse_trip_mpz
(
    SLIP_sparse **A_handle,     // matrix to construct
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    mpz_t *x,           // Set of values in full precision integer
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x, I, and J vectors)
)
{
    SLIP_info ok;
    if (!I || !J || !A_handle || !x || n <= 0 || nz <= 0)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_sparse *A = slip_create_sparse ( ) ;
    if (A == NULL)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    SLIP_CHECK(SLIP_mpq_set_ui(A->scale, 1, 1));

    SLIP_CHECK(slip_trip_to_mat(A, I, J, x, n, nz));

    (*A_handle) = A ;
    return SLIP_OK;
}
