//------------------------------------------------------------------------------
// SLIP_LU/slip_back_sub: sparse REF backward substitution (x = U\x)
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse REF backward substitution, solving
 * the sysem Ux = b. Note that prior to this, we expect x to be multiplied by
 * the determinant of A.
 *
 * U is a sparse mpz matrix, and bx is a dense mpz matrix.  The diagonal entry
 * of U must appear as the last entry in each column.
 *
 * The input argument bx contains b on input, and it is overwritten on output
 * by the solution x.
 */

#include "SLIP_LU_internal.h"

SLIP_info slip_back_sub  // performs sparse REF backward substitution
(
    const SLIP_sparse *U,   // input upper triangular matrix
    mpz_t **bx,             // right hand side matrix of size n*numRHS
    int64_t numRHS          // number of columns in bx
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (U,  SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (bx, SLIP_DENSE, SLIP_MPZ) ;

    //--------------------------------------------------------------------------

    int sgn;
    mpz_t *Ux = U->x;
    int64_t *Ui = U->i;
    int64_t *Up = U->p;

    for (int64_t k = 0; k < numRHS; k++)
    {
        // Start at bx[n]
        for (int64_t j = U->n-1; j >= 0; j--)
        {
            // If bx[j] is zero skip this iteration
                SLIP_CHECK(SLIP_mpz_sgn(&sgn, bx[j][k]));
            if (sgn == 0) {continue;}

            // Obtain bx[j]
            SLIP_CHECK(SLIP_mpz_divexact(bx[j][k], bx[j][k],Ux[Up[j+1]-1]));
            for (int64_t i = Up[j]; i < Up[j+1]-1; i++)
            {
                SLIP_CHECK(SLIP_mpz_sgn(&sgn, Ux[i]));
                if (sgn == 0) {continue;}
                // bx[i] = bx[i] - Ux[i]*bx[j]
                SLIP_CHECK(SLIP_mpz_submul(bx[Ui[i]][k], Ux[i], bx[j][k]));
            }
        }
    }

    return (SLIP_OK) ;
}

