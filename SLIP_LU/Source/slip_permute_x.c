//------------------------------------------------------------------------------
// SLIP_LU/slip_permute_x: permute x, as x = Q*x
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function permutes x to get it back in its original form.
 * That is, x = Q*x.
 */

#define SLIP_FREE_ALL                 \
    SLIP_matrix_free(&x2, NULL);

#include "SLIP_LU_internal.h"

SLIP_info slip_permute_x
(
    SLIP_matrix *x,       // Solution vector
    SLIP_LU_analysis *S,  // symbolic analysis with the column ordering Q
    SLIP_options* option  // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (x, SLIP_DENSE, SLIP_MPQ) ;

    if (!x || !S || !S->q) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------
    // x2 (q) = x
    //--------------------------------------------------------------------------

    int64_t *q = S->q, n = x->m, numRHS = x->n;     // column permutation
    // Declare temp x
    SLIP_matrix *x2;
    SLIP_matrix_allocate(&x2, SLIP_DENSE, SLIP_MPQ, x->m, x->n, x->nzmax,
                         false, true, option);
    if (!x2) {return SLIP_OUT_OF_MEMORY;}
    // Set x2 = Q*x
    for (int64_t i = 0; i < n; i++)
    {
        for (int64_t j = 0; j < numRHS; j++)
        {
            SLIP_CHECK(SLIP_mpq_set(SLIP_2D(x2, q[i], j, mpq),
                                    SLIP_2D(x, i, j, mpq)));
        }
    }

    //--------------------------------------------------------------------------
    // x = x2
    //--------------------------------------------------------------------------

    for (int64_t i = 0; i < n; i++)
    {
        for (int64_t j = 0; j < numRHS; j++)
        {
            SLIP_CHECK(SLIP_mpq_set(SLIP_2D(x, i, j, mpq),
                                    SLIP_2D(x2, i, j, mpq)));
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE_ALL;
    return SLIP_OK;
}

