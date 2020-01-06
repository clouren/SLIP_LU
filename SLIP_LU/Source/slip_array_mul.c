//------------------------------------------------------------------------------
// SLIP_LU/slip_array_mul: multiplies a vector by a scalar
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* 
 * Purpose: This function multiplies vector x by the determinant of matrix. 
 * 
 * On output the contents of the x vector is modified 
 */

#include "SLIP_LU_internal.h"

SLIP_info slip_array_mul // multiplies vector x by the determinant of matrix
(
    mpz_t** x,      // matrix to be multiplied
    mpz_t det,      // given determinant of matrix
    int32_t n,      // size of x
    int32_t numRHS  // number of RHS vectors
)
{
    SLIP_info ok = SLIP_OK;
    if (!x) {return SLIP_INCORRECT_INPUT;}
    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t k = 0; k < numRHS; k++)
        {
            // x[i][k] = x[i][k] * det
            SLIP_CHECK(SLIP_mpz_mul(x[i][k], x[i][k], det));
        }
    }
    return ok;
}

