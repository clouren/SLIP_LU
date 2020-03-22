//------------------------------------------------------------------------------
// SLIP_LU/SLIP_check_solution: check solution to Ax=b
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Check the solution of the linear system by performing a quick
 * rational arithmetic A*x = b.
 *
 * WARNING: This function WILL produce incorrect results if called with the
 * final scaled x vector! This function assumes that A and b are in their
 * scaled integral form! Thus x is the solution to the scaled system Ax=b. Once
 * x has been finalized, this function will report an incorrect result. Please
 * refer to Demo/SLIPLU.c for the proper usage of this function.
 */

#define SLIP_FREE_ALL                       \
    SLIP_MPQ_CLEAR(temp);                   \
    SLIP_delete_mpq_mat(&b2, n, numRHS);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_check_solution
(
    SLIP_sparse *A,          // input matrix
    mpq_t** x,               // solution vector
    SLIP_dense *b            // right hand side
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (A, SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (x, SLIP_DENSE, SLIP_MPQ) ;
    SLIP_REQUIRE (b, SLIP_DENSE, SLIP_MPZ) ;

    if (!b->x || !A->p || !A->i || !A->x)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare vars
    //--------------------------------------------------------------------------

    int64_t p, j, i ;
    int64_t n = A->n, numRHS = b->n;
    mpz_t** bx = b->x;
    mpq_t** b2 = NULL;
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(SLIP_mpq_init(temp));

    // b2 stores the solution of A*x
    b2 = SLIP_create_mpq_mat(n, numRHS);
    if (!b2)
    {
        SLIP_FREE_ALL;
        return SLIP_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // perform SLIP_mpq_addmul in loops
    //--------------------------------------------------------------------------

    for (j = 0; j < numRHS; j++)
    {
        for (i = 0; i < n; i++)
        {
            for (p = A->p[i]; p < A->p[i + 1]; p++)
            {
                // temp = A[p][i]
                SLIP_CHECK(SLIP_mpq_set_z(temp, A->x[p]));

                // temp = temp*x[i]
                SLIP_CHECK(SLIP_mpq_mul(temp, temp, x[i][j]));

                // b2[p] = b2[p]-temp
                SLIP_CHECK(SLIP_mpq_add(b2[A->i[p]][j], b2[A->i[p]][j],temp));
            }
        }
    }

    //--------------------------------------------------------------------------
    // check if b==b2
    //--------------------------------------------------------------------------

    for (j = 0; j < numRHS; j++)
    {
        for (i = 0; i < n; i++)
        {
            // z = b[i] (correct b)
            SLIP_CHECK(SLIP_mpq_set_z(temp, bx[i][j]));

            // set check false if b!=b2
            int r ;
            SLIP_CHECK(SLIP_mpq_equal(&r, temp, b2[i][j]));
            if (r == 0)
            {
                SLIP_FREE_ALL;
                return SLIP_INCORRECT;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    SLIP_FREE_ALL;
    return SLIP_OK;
}

