# include "SLIP_LU_internal.h"


/* 
 * Purpose: This function multiplies vector x by the determinant of matrix. 
 * 
 * On output the contents of the x vector is modified 
 */
#define SLIP_FREE_WORKSPACE

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
            SLIP_CHECK(slip_mpz_mul(x[i][k], x[i][k], det));
        }
    }
    return ok;
}
#undef SLIP_FREE_WORKSPACE
