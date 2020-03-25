//------------------------------------------------------------------------------
// SLIP_LU/slip_scale_x: scale the matrix x
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function scales the x matrix if necessary */

#include "SLIP_LU_internal.h"

SLIP_info slip_scale_x
(
    SLIP_matrix *x,         // Solution matrix
    SLIP_matrix *A,         // matrix A
    SLIP_matrix *b,         // right hand side
    SLIP_options* option    // command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (A, SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (b, SLIP_DENSE, SLIP_MPZ) ;
    SLIP_REQUIRE (x, SLIP_DENSE, SLIP_MPQ) ;
    if (!option)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------

    int r, s ;
    int64_t n, numRHS;
    n = A->m;
    numRHS = b->n;

    // Determine if A's scaling factor is 1
    SLIP_CHECK(SLIP_mpq_cmp_ui(&r, A->scale, 1, 1));
    SLIP_CHECK(SLIP_mpq_cmp_ui(&s, A->scale, 0, 1));
    if (r != 0 && s != 0)
    {
        for (int64_t i = 0; i < n; i++)
        {
            for (int64_t j = 0; j < numRHS; j++)
            {
                SLIP_CHECK(SLIP_mpq_mul(SLIP_2D(x, i, j, mpq),
                                        SLIP_2D(x, i, j, mpq),
                                        A->scale));
            }
        }
    }

    // Determine if b's scaling factor is 1
    SLIP_CHECK(SLIP_mpq_cmp_ui(&r, b->scale, 1, 1));
    SLIP_CHECK(SLIP_mpq_cmp_ui(&s, b->scale, 0, 1));
    if (r != 0 && s != 0)
    {
        for (int64_t i = 0; i < n; i++)
        {
            for (int64_t j = 0; j < numRHS; j++)
            {
                SLIP_CHECK(SLIP_mpq_div(SLIP_2D(x, i, j, mpq), 
                                        SLIP_2D(x, i, j, mpq),
                                        b->scale));
            }
        }
    }
    return SLIP_OK;
}
