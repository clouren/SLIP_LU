//------------------------------------------------------------------------------
// SLIP_LU/SLIP_LU_solve: exact solution of Ax=b
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system LD^(-1)U x = b. It
 * essnetially serves as a wrapper for all forward and backward substitution
 * routines. This function always returns the solution matrix x as a rational
 * matrix. If a user desires to have double or mpfr output, they must create
 * a matrix copy.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Unitialized on input.
 *           on output, contains the exact rational solution of the system
 *
 * b:        Set of RHS vectors
 *
 * A:        Input matrix. Unmodified on input/output
 *
 * L:        Lower triangular matrix. Unmodified on input/output
 *
 * U:        Upper triangular matrix. Unmodified on input/output
 *
 * rhos:     dense mpz_t matrix of pivots. Contains the sequence of pivots
 *           encountered during factorization and is used for forward/back
 *           substitution. Unmodified on input/output.
 *
 * S:        symbolic analysis struct, used for vector permutation
 *
 * pinv:     inverse row permutation vector, used to permute the b vectors.
 *           unmodified on input/output.
 *
 * option:   command options
 */

#define SLIP_FREE_WORK                        \
    SLIP_matrix_free(&b2, NULL);              \
    SLIP_MPQ_CLEAR(scale);                    \

#define SLIP_FREE_ALL                         \
    SLIP_FREE_WORK                            \
    SLIP_matrix_free(&x, NULL);               \

#include "slip_LU_internal.h"

SLIP_info SLIP_LU_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SLIP_matrix **x_handle,  // rational solution to the system
    // input:
    const SLIP_matrix *b,   // right hand side vector
    const SLIP_matrix *A,   // Input matrix
    const SLIP_matrix *L,   // lower triangular matrix
    const SLIP_matrix *U,   // upper triangular matrix
    const SLIP_matrix *rhos,// sequence of pivots
    const SLIP_LU_analysis *S,// symbolic analysis struct
    const int64_t *pinv,    // row permutation
    const SLIP_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------


    SLIP_info info ;
    SLIP_REQUIRE (b,    SLIP_DENSE, SLIP_MPZ) ;
    SLIP_REQUIRE (A,    SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (L,    SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (U,    SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (rhos, SLIP_DENSE, SLIP_MPZ) ;

    if (!x_handle || !S || !pinv || L->m != A->m || L->n != U->m ||
        U->n != A->n || A->n != A->m || A->m != b->m )
    {
        return SLIP_INCORRECT_INPUT;
    }
    *x_handle = NULL;


    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    int64_t i, k, n = L->n, numRHS = b->n;
    mpq_t scale;
    SLIP_MPQ_SET_NULL(scale);

    SLIP_matrix *x = NULL;   // solution
    SLIP_matrix *b2 = NULL;  // permuted b
    SLIP_CHECK (SLIP_matrix_allocate(&x, SLIP_DENSE, SLIP_MPQ, n, numRHS,
        n*numRHS, false, true, option));

    // not need to perform deep copy, simply create an new dense matrix
    SLIP_CHECK (SLIP_matrix_allocate(&b2, SLIP_DENSE, SLIP_MPZ, n, numRHS,
        n*numRHS, false, false, option));

    //--------------------------------------------------------------------------
    // Set workspace b
    //--------------------------------------------------------------------------

    for (k = 0; k < numRHS; k++)
    {
        for (i = 0; i < n; i++)
        {
            SLIP_CHECK(SLIP_mpz_set(SLIP_2D(b2, pinv[i], k, mpz),
                                    SLIP_2D(b,  i,       k, mpz)));
        }
    }

    //--------------------------------------------------------------------------
    // b2 = L\b2, via forward substitution
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_forward_sub(L, b2, (SLIP_matrix*) rhos));

    //--------------------------------------------------------------------------
    // b2 = b2 * det, where det=rhos[n-1]
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_matrix_mul(b2, rhos->x.mpz[n-1]));

    //--------------------------------------------------------------------------
    // b2 = U\b2, via back substitution
    //--------------------------------------------------------------------------
    SLIP_CHECK(slip_back_sub(U, b2));

    //--------------------------------------------------------------------------
    // x = b2/det, where det=rhos[n-1]
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_matrix_div(x, b2, rhos->x.mpz[n-1]));

    //--------------------------------------------------------------------------
    // Permute the solution vectors
    //--------------------------------------------------------------------------

    // TODO Please refer to the comment in the below function. Once a solution
    // is decided update the rest of this file.
    SLIP_CHECK(slip_permute_x(x, (SLIP_LU_analysis *) S, option));

    //--------------------------------------------------------------------------
    // Check the solution if desired (debugging only)
    //--------------------------------------------------------------------------

    if (option->check)
    {
        SLIP_CHECK(slip_check_solution(A, x, b, option));
    }

    //--------------------------------------------------------------------------
    // Scale the solution if necessary.
    //--------------------------------------------------------------------------

    SLIP_CHECK(SLIP_mpq_init(scale));

    // set the scaling factor scale = A->scale / b->scale
    SLIP_CHECK( SLIP_mpq_div(scale, A->scale, b->scale));

    // Determine if the scaling factor is 1
    int r;
    SLIP_CHECK(SLIP_mpq_cmp_ui(&r, scale, 1, 1));
    int64_t nz = x->m * x->n;
    if (r != 0 )
    {
        for (i = 0; i < nz; i++)
        {
            SLIP_CHECK(SLIP_mpq_mul(x->x.mpq[i], x->x.mpq[i], scale));
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE_WORK;
    (*x_handle) = x;
    return SLIP_OK;
}

