//------------------------------------------------------------------------------
// SLIP_LU/slip_back_sub: sparse REF backward substitution (x = U\x)
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse REF backward substitution. In essense
 * it solves the sysem Ux = x. Note that prior to this, we expect x to be
 * multiplied by the determinant of A.
 *
 * The input argument x is modified on output
 *
 */

#include "SLIP_LU_internal.h"

SLIP_info slip_back_sub  // performs sparse REF backward substitution
(
    SLIP_sparse *U,   // input upper triangular matrix
    mpz_t **bx,       // right hand side matrix of size n*numRHS
    int32_t numRHS    // number of columns in bx
)
{
    int32_t sgn;
    mpz_t *Ux = U->x;
    int32_t *Ui = U->i;
    int32_t *Up = U->p;
    SLIP_info ok;
    for (int32_t k = 0; k < numRHS; k++)
    {
        // Start at bx[n]
        for (int32_t j = U->n-1; j >= 0; j--)
        {
            // If bx[j] is zero skip this iteration
                SLIP_CHECK(SLIP_mpz_sgn(&sgn, bx[j][k]));
            if (sgn == 0) {continue;}

            // Obtain bx[j]
            SLIP_CHECK(SLIP_mpz_divexact(bx[j][k], bx[j][k],Ux[Up[j+1]-1]));
            for (int32_t i = Up[j]; i < Up[j+1]-1; i++)
            {
                SLIP_CHECK(SLIP_mpz_sgn(&sgn, Ux[i]));
                if (sgn == 0) {continue;}
                // bx[i] = bx[i] - Ux[i]*bx[j]
                SLIP_CHECK(SLIP_mpz_submul(bx[Ui[i]][k], Ux[i], bx[j][k]));
            }
        }
    }
    return SLIP_OK;
}

