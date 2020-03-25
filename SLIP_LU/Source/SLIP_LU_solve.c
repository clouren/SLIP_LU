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
 * X_handle: A pointer to the solution vectors. Unitialized on input.
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

#define SLIP_FREE_ALL                        \
    SLIP_matrix_free(&b2, NULL);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_LU_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SLIP_matrix **X_handle,  // rational solution to the system
    // input:
    SLIP_matrix *b,         // right hand side vector
    const SLIP_matrix *A,   // Input matrix
    const SLIP_matrix *L,   // lower triangular matrix
    const SLIP_matrix *U,   // upper triangular matrix
    const SLIP_matrix *rhos,// sequence of pivots
    SLIP_LU_analysis *S,    // symbolic analysis struct
    const int64_t *pinv,    // row permutation
    SLIP_options* option    // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (A,    SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (L,    SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (U,    SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (rhos, SLIP_DENSE, SLIP_MPZ) ;
    SLIP_REQUIRE (b,    SLIP_DENSE, SLIP_MPZ) ;
    if (!option)
        return SLIP_INCORRECT_INPUT;

    SLIP_matrix *x = NULL;
    SLIP_matrix_allocate(&x, SLIP_DENSE, SLIP_MPQ, b->m, b->n, b->m*b->n, false, true, option);
    
    if (!x || !b || !rhos || !pinv || !L || !U 
        || !L->p || !L->i || !U->p || !U->i)
    {
        return SLIP_INCORRECT_INPUT;
    }


    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    int64_t i, k, n = L->n;
    
    //--------------------------------------------------------------------------
    // create the permuted b
    //--------------------------------------------------------------------------

    SLIP_matrix *b2;
    
    SLIP_matrix_copy(&b2, SLIP_DENSE, SLIP_MPZ, b, option);
    
    if (!b2)
    {
        return SLIP_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // Set workspace b
    //--------------------------------------------------------------------------

    for (k = 0; k < b->n; k++)
    {
        for (i = 0; i < b->m; i++)
        {
            SLIP_CHECK(SLIP_mpz_set(SLIP_2D(b2, pinv[i], k, mpz), 
                                    SLIP_2D(b, i, k, mpz)));
        }
    }

    //--------------------------------------------------------------------------
    // b2 = L2\b2, via forward substitution
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_forward_sub(L, b2, (SLIP_matrix*) rhos, option));

    //--------------------------------------------------------------------------
    // b2 = b2 * det, where det=rhos[n-1]
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_array_mul(b2, rhos->x.mpz[n-1], option));

    //--------------------------------------------------------------------------
    // b2 = U\b2, via back substitution
    //--------------------------------------------------------------------------
    SLIP_CHECK(slip_back_sub(U, b2, option));

    //--------------------------------------------------------------------------
    // x = b2/det, where det=rhos[n-1]
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_array_div(x, b2, rhos->x.mpz[n-1], option));

    //--------------------------------------------------------------------------
    // Permute the solution vectors
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_permute_x(x, S, option));
    
    //--------------------------------------------------------------------------
    // Check the solution if desired (debugging only)
    //--------------------------------------------------------------------------

    if (option->check)
    {
        SLIP_info checker = slip_check_solution( (SLIP_matrix*) A, x, b, option);
        if (checker == SLIP_OK)
        {
            printf ("Solution is verified to be exact.\n") ;
        }
        else if (checker == SLIP_INCORRECT)
        {
            // This can never happen.
            printf ("ERROR! Solution is wrong. This is a bug; please"
                "contact the authors of SLIP LU.\n") ;
            abort ( ) ;
        }
        else
        {
            // Out of memory or bad input.
            SLIP_FREE_ALL;
            return checker;
        }
    }
            
    //--------------------------------------------------------------------------
    // scale solution
    //--------------------------------------------------------------------------

    SLIP_CHECK(slip_scale_x(x, (SLIP_matrix*) A, b, option));
    
    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*X_handle) = x;
    SLIP_FREE_ALL;
    return SLIP_OK;
}

