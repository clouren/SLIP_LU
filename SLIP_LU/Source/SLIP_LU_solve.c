//------------------------------------------------------------------------------
// SLIP_LU/SLIP_LU_solve: exact solution of Ax=b
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system LD^(-1)U x = b.*/

#define SLIP_FREE_WORKSPACE                        \
    SLIP_delete_mpz_mat(&b2, n, numRHS);

# include "SLIP_LU_internal.h"

SLIP_info SLIP_LU_solve  //solves the linear system LD^(-1)U x = b
(
    mpq_t **x,           // rational solution to the system
    SLIP_dense *b,       // right hand side vector
    mpz_t *rhos,         // sequence of pivots
    SLIP_sparse *L,      // lower triangular matrix
    SLIP_sparse *U,      // upper triangular matrix
    int32_t *pinv        // row permutation
)
{
    if (!x || !b || !rhos || !pinv || !L || !U || !b->x
        || !L->p || !L->i || !L->x || !U->p || !U->i || !U->x)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t i, k, n = L->n, numRHS = b->n;
    mpz_t **bx = b->x;
    SLIP_info ok;

    // Permuted b
    mpz_t** b2 = SLIP_create_mpz_mat(n, numRHS);
    if (!b2)
    {
        return SLIP_OUT_OF_MEMORY;
    }

    // Set workspace b
    for (k = 0; k < numRHS; k++)
    {
        for (i = 0; i < n; i++)
        {
            SLIP_CHECK(SLIP_mpz_set(b2[pinv[i]][k], bx[i][k]));
        }
    }
    // L*b2 = b2
    SLIP_CHECK(slip_forward_sub(L, b2, rhos, numRHS));
    // b2 = b2 * det, where det=rhos[n-1]
    SLIP_CHECK(slip_array_mul(b2, rhos[n-1], n, numRHS));
    // U b2 = b2
    SLIP_CHECK(slip_back_sub(U, b2, numRHS));
    // x = b2/det, where det=rhos[n-1]
    SLIP_CHECK(slip_array_div(x, b2, rhos[n-1], n, numRHS));
    // Free memory
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

