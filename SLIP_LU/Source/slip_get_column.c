//------------------------------------------------------------------------------
// SLIP_LU/slip_get_column: extract a column from a matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function obtains column k from matrix A and stores it in the 
 * dense vector x
 * 
 * On exit, x either contains the kth column of A or is NULL
 */

#define SLIP_FREE_WORKSPACE                \
    SLIP_delete_mpz_array(&x, A->n);

# include "SLIP_LU_internal.h"

SLIP_info slip_get_column //extract k-th column from A, i.e., x=A(:,k)
(
    mpz_t* x,       // A(:,k)
    SLIP_sparse* A, // input matrix
    int32_t k       // column to extract
)
{
    if (!A || !x || !A->x || !A->p || !A->i) {return SLIP_INCORRECT_INPUT;}
    SLIP_info ok = SLIP_OK;
    // Iterating accross the nonzeros in column k
    for (int32_t i = A->p[k]; i < A->p[k + 1]; i++)
    {
        // Value of the ith nonzero
        SLIP_CHECK(SLIP_mpz_set(x[A->i[i]], A->x[i]));
    }
    return SLIP_OK;
}

