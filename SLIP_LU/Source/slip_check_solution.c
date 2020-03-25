//------------------------------------------------------------------------------
// SLIP_LU/slip_check_solution: check solution to Ax=b
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
 * scaled integral form! Thus x is the solution to the scaled system Ax=b. Please
 * note the usage in SLIP_LU_solve
 */

#define SLIP_FREE_ALL                       \
    SLIP_MPQ_CLEAR(temp);                   \
    SLIP_matrix_free(&b2, NULL);

#include "SLIP_LU_internal.h"

SLIP_info slip_check_solution
(
    SLIP_matrix *A,         // Input matrix
    SLIP_matrix *x,         // Solution vectors
    SLIP_matrix *b,         // Right hand side vectors
    SLIP_options* option    // Command options, currently unused
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (A, SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (x, SLIP_DENSE, SLIP_MPQ) ;
    SLIP_REQUIRE (b, SLIP_DENSE, SLIP_MPZ) ;

    if ( !A->p || !A->i || !option)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare vars
    //--------------------------------------------------------------------------

    int64_t p, j, i ;
    SLIP_matrix *b2 = NULL;
     SLIP_matrix_allocate(&b2, SLIP_DENSE, SLIP_MPQ, b->m, b->n, b->nzmax, 
                         false, true, option);

    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(SLIP_mpq_init(temp));

     // b2 stores the solution of A*x
    if (!b2)
    {
        SLIP_FREE_ALL;
        return SLIP_OUT_OF_MEMORY;
    }

     //--------------------------------------------------------------------------
    // perform SLIP_mpq_addmul in loops
    //--------------------------------------------------------------------------
    
    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            for (p = A->p[i]; p < A->p[i + 1]; p++)
            {
                // temp = A[p][i]
                SLIP_CHECK(SLIP_mpq_set_z(temp, A->x.mpz[p]));

                // temp = temp*x[i]
                SLIP_CHECK(SLIP_mpq_mul(temp, temp, 
                                        SLIP_2D(x, i, j, mpq)));

                // b2[p] = b2[p]-temp
                SLIP_CHECK(SLIP_mpq_add(SLIP_2D(b2, A->i[p], j, mpq),
                                        SLIP_2D(b2, A->i[p], j, mpq),temp));
            }
        }
    }

    //--------------------------------------------------------------------------
    // check if b==b2
    //--------------------------------------------------------------------------

    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            // z = b[i] (correct b)
            SLIP_CHECK(SLIP_mpq_set_z(temp, SLIP_2D(b, i, j, mpz)));

            // set check false if b!=b2
            int r ;
            SLIP_CHECK(SLIP_mpq_equal(&r, temp, SLIP_2D(b2, i, j, mpq)));
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

