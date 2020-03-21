//------------------------------------------------------------------------------
// SLIP_LU/SLIP_permute_x: permute x, as x = Q*x
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function permutes x to get it back in its original form.
 * That is, x = Q*x.
 */

#define SLIP_FREE_ALL                 \
    SLIP_delete_mpq_mat(&x2,n,numRHS);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_permute_x
(
    mpq_t **x,            // Solution vector
    int32_t n,            // Size of solution vector
    int32_t numRHS,       // number of RHS vectors
    SLIP_LU_analysis *S   // symbolic analysis with the column ordering Q
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: change this to create a new SLIP_dense matrix on output?
    // It would be faster, since the copy x = x2 could be skipped.

    SLIP_info info ;
    SLIP_REQUIRE (x, SLIP_DENSE, SLIP_MPQ) ;

    if (!x || !S || !S->q) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------
    // x2 (q) = x
    //--------------------------------------------------------------------------

    int32_t *q = S->q;     // column permutation
    // Declare temp x
    mpq_t** x2 = SLIP_create_mpq_mat(n, numRHS);
    if (!x2) {return SLIP_OUT_OF_MEMORY;}
    // Set x2 = Q*x
    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t j = 0; j < numRHS; j++)
        {
            SLIP_CHECK(SLIP_mpq_set(x2[q[i]][j], x[i][j]));
        }
    }

    //--------------------------------------------------------------------------
    // x = x2
    //--------------------------------------------------------------------------

    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t j = 0; j < numRHS; j++)
        {
            SLIP_CHECK(SLIP_mpq_set(x[i][j], x2[i][j]));
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE_ALL;
    return SLIP_OK;
}

