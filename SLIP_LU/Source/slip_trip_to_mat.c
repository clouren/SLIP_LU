//------------------------------------------------------------------------------
// SLIP_LU/slip_trip_to_mat: convert triplet to sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/*
 * Purpose: This function converts triplet matrix into compressed column
 * matrix A
 *
 * The values of the triplet matrix must be stored as mpz_t values
 */

#define SLIP_FREE_WORKSPACE              \
    SLIP_FREE(w);

#include "SLIP_LU_internal.h"

SLIP_info slip_trip_to_mat
(
    SLIP_sparse *A,     //matrix stored in ccf that will take B
    int32_t *I,         // Row indices.
    int32_t *J,         // Column indices
    mpz_t *x,           // Values in the matrix
    int32_t n,          // Dimension of the matrix
    int32_t nz          // Number of nonzeros in the matrix
)
{
    // inputs have been validated in SLIP_build_sparse_trip_*.c
    int32_t k, p;
    SLIP_info ok;
    int32_t* w = (int32_t*) SLIP_calloc(n, sizeof(int32_t));
    if (!w)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    SLIP_CHECK(slip_sparse_alloc(A, n, n, nz));

    // Column pointers
    for (k = 0; k < nz; k++)
    {
        w[J[k]]++;
    }
    // Column sums for A, w = A->p[0..n-1]
    slip_cumsum(A->p, w, A->n);
    for (k = 0; k < nz; k++)
    {
        p = w[J[k]]++;
        // Place values of i
        A->i[p] = I[k];
        // Place values of x
        if (A->i[p] < 0)
        {
            SLIP_FREE_WORKSPACE;
            return SLIP_INCORRECT_INPUT;
        }
        SLIP_CHECK(SLIP_mpz_set(A->x[p], x[k]));
    }
    // Number of nonzeros in A
    A->nz = A->nzmax;
    // Delete w
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

