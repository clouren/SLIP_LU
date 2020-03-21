//------------------------------------------------------------------------------
// SLIP_LU/slip_dense_alloc: allocate a dense matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function allocates a SLIP dense matrix of size n*m.  This
 * function manually initializes each entry in A->x therefore they are
 * immediately ready to be operated on.
 */

#include "SLIP_LU_internal.h"

SLIP_info slip_dense_alloc
(
    SLIP_dense* A, // dense matrix data structure to be allocated
    int64_t m,     // number of rows
    int64_t n      // number of columns
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_allocate (&A, SLIP_DENSE, SLIP_MPZ, ...) ;

    if (n <= 0 || m <= 0)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------

    A->m = m;                                    // Rows of the matrix
    A->n = n;                                    // Columns of the matrix
    A->x = SLIP_create_mpz_mat(m, n);            // Initialize the mpz mat x
    if (!A->x)
    {
        return SLIP_OUT_OF_MEMORY;
    }
    return SLIP_OK;
}
