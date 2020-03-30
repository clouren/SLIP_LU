//------------------------------------------------------------------------------
// SLIP_LU/slip_array_mul: multiplies a vector by a scalar
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function multiplies the dense matrix x by the determinant of matrix.
 * 
 * This function requires that the matrix x is mpz_t and dense
 *
 * On output the contents of the x vector is modified.
 */

#include "SLIP_LU_internal.h"

// TODO rename slip_matrix_mul?
SLIP_info slip_array_mul // multiplies vector x by the determinant of matrix
(
    SLIP_matrix *x,         // matrix to be multiplied
    const mpz_t det,        // given determinant of matrix
    SLIP_options* option    // Command options, currently unused
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (x, SLIP_DENSE, SLIP_MPZ) ;
    if (!option)
    {
        return SLIP_INCORRECT_INPUT;// TODO create default option?
    }
    //--------------------------------------------------------------------------
    // x = x * det
    //--------------------------------------------------------------------------

    for (int64_t i = 0; i < x->n * x->m; i++)
    {
        // x[i] = x[i]*det
        SLIP_CHECK( SLIP_mpz_mul( x->x.mpz[i], x->x.mpz[i], det));
    }

    return (SLIP_OK) ;
}

