//------------------------------------------------------------------------------
// SLIP_LU/slip_get_column: extract a column from a matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function obtains column k from matrix A and scatters it into
 * the dense vector x.
 *
 * On exit, x either contains the kth column of A or is NULL
 */

#define SLIP_FREE_ALL                \
    SLIP_delete_mpz_array(&x, A->n);

#include "SLIP_LU_internal.h"

SLIP_info slip_get_column //extract k-th column from A, i.e., x=A(:,k)
(
    mpz_t* x,       // A(:,k)
    SLIP_sparse* A, // input matrix
    int64_t k       // column to extract
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: this function is only used in slip_ref_triangular_solve.
    // Just delete this function and move the functionality there?

    SLIP_info info ;
    ASSERT (A != NULL && A->kind == SLIP_CSC && A->type == SLIP_MPZ) ;
    ASSERT (x != NULL) ;        // TODO make x a SLIP_matrix?

    if (!A || !x || !A->x || !A->p || !A->i) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------

    // Iterating accross the nonzeros in column k
    for (int64_t i = A->p[k]; i < A->p[k + 1]; i++)
    {
        // Value of the ith nonzero
        SLIP_CHECK(SLIP_mpz_set(x[A->i[i]], A->x[i]));
    }
    return SLIP_OK;
}

