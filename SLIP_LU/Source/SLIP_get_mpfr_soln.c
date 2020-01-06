//------------------------------------------------------------------------------
// SLIP_LU/SLIP_get_mpfr_soln: convert mpq solution to mpfr
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

/* Purpose: Convert the output mpq_t** solution vector obtained from
 * SLIP_LU_solve and SLIP_permute_x from mpq_t** to mpfr_t**.
 * x_mpfr has to be initialized before passed in.
 */

SLIP_info SLIP_get_mpfr_soln
(
    mpfr_t **x_mpfr,      // mpfr solution of size n*numRHS to Ax = b
    mpq_t  **x_mpq,       // mpq solution of size n*numRHS to Ax = b.
    int32_t n,            // Dimension of A, number of rows of x 
    int32_t numRHS,        // Number of right hand side vectors
    SLIP_options* option
)
{
    if (x_mpfr  == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t j = 0; j < numRHS; j++)
        {
            SLIP_info ok = SLIP_mpfr_set_q(x_mpfr[i][j], x_mpq[i][j],
                option->SLIP_MPFR_ROUND);
            if (ok != SLIP_OK)         {  return ok;  }
        }
    }
    return SLIP_OK;
}

