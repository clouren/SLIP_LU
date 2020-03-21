//------------------------------------------------------------------------------
// SLIP_LU/slip_mpz_populate_mat: construct a sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function populates the SLIP_sparse A by the CSC stored vectors
 * I, p, and x.
 */

#include "SLIP_LU_internal.h"

SLIP_info slip_mpz_populate_mat
(
    SLIP_sparse* A,   // matrix to be populated
    int64_t* I,       // row indices
    int64_t* p,       // column pointers
    mpz_t* x,         // set of values in A
    int64_t n,        // size of the matrix A
    int64_t nz        // number of nonzeros in the matrix A
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: delete this.  It becomes part of SLIP_matrix_copy.

    // inputs have been validated in SLIP_build_sparse_csc_*.c
    SLIP_info info ;

    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_sparse_alloc(A, n, n, nz));// Allocate the matrix A
    A->nz = nz;
    for (int64_t k = 0; k <= n; k++)              // Set A->p
    {
        A->p[k] = p[k];
    }
    for (int64_t k = 0; k < nz; k++)              // Set A->i and A->x
    {
        A->i[k] = I[k];
        if (A->i[k] < 0)
        {
            return SLIP_INCORRECT_INPUT;
        }
        SLIP_CHECK(SLIP_mpz_set(A->x[k],x[k]));
    }

    return (SLIP_OK) ;
}

