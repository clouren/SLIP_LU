//------------------------------------------------------------------------------
// SLIP_LU/SLIP_solve_mpfr: solve Ax=b, returning solution as mpfr matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------


/* This code utilizes the SLIP LU factorization. Soln is output as mpfr_t.
 *
 * Input/Output arguments:
 *
 * x_mpfr:  A mpfr_t** array which is unitialized on input. On output, contains
 *          the solution to the linear system Ax=b in user specified precision.
 *
 * A:       User's input matrix. It must be populated prior to calling this
 *          function.
 *
 * S:       Symbolic analysis struct. It is important that this data structure
 *          is both allocated and populated prior to calling this function.
 *          That is, this function MUST be preceded by a call to the function
 *          SLIP_LU_analyze.
 *
 * b:       Collection of right hand side vectors. Must be populated prior to
 *          factorization.
 *
 * option:  Struct containing various command parameters for the factorization.
 */

# define SLIP_FREE_ALL                \
    SLIP_delete_sparse(&L);                 \
    SLIP_delete_sparse(&U);                 \
    SLIP_FREE(pinv);                        \
    SLIP_delete_mpz_array(&rhos, n);        \
    SLIP_delete_mpq_mat(&x_mpq, n, numRHS);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_solve_mpfr
(
    mpfr_t **x_mpfr,        // Solution vector stored as an mpfr_t array
    SLIP_sparse *A,         // Compressed column form full precision matrix A
    SLIP_LU_analysis *S,    // Column ordering
    SLIP_dense *b,          // Right hand side vectrors
    SLIP_options *option    // Control parameters
)
{

    //-------------------------------------------------------------------------
    // check inputs
    //-------------------------------------------------------------------------

    SLIP_info info ;
    if (!x_mpfr || !A || !A->p || !A->i || !A->x ||
        !S || !S->q || !b || !b->x || !option)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_sparse *L = NULL ;
    SLIP_sparse *U = NULL ;
    int32_t *pinv = NULL ;
    mpz_t *rhos = NULL ;

    //--------------------------------------------------------------------------
    // Declare memory
    //--------------------------------------------------------------------------

    int32_t n = A->n, numRHS = b->n;
    mpq_t **x_mpq = SLIP_create_mpq_mat(n, numRHS);
    if (!x_mpq)
    {
        SLIP_FREE_ALL;
        return SLIP_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // LU factorization
    //--------------------------------------------------------------------------

    SLIP_CHECK(SLIP_LU_factorize(&L, &U, &rhos, &pinv, A, S, option));

    //--------------------------------------------------------------------------
    // FB substituion
    //--------------------------------------------------------------------------

    SLIP_CHECK (SLIP_LU_solve(x_mpq, b,
        (const SLIP_sparse *) L,
        (const SLIP_sparse *) U,
        (const mpz_t *) rhos,
        (const int32_t *) pinv)) ;

    SLIP_CHECK(SLIP_permute_x(x_mpq, n, numRHS, S));

#if 0
    // Check solution (debugging only)
    SLIP_CHECK(SLIP_check_solution(A, x_mpq, b));
#endif

    //--------------------------------------------------------------------------
    // Scale solution
    //--------------------------------------------------------------------------

    SLIP_CHECK(SLIP_scale_x(x_mpq, A, b));

    //--------------------------------------------------------------------------
    // Output and free memory
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_get_mpfr_soln(x_mpfr, x_mpq, n, numRHS, option));

    SLIP_FREE_ALL;
    return (SLIP_OK) ;
}

