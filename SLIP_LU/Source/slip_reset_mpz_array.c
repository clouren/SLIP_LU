//------------------------------------------------------------------------------
// SLIP_LU/slip_reset_mpz_array: clear the used entries of an mpz workspace
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function resets an mpz array of size n with the nonzero
 * pattern given. This is more efficient than iterating accross all nonzeros in
 * vector x.
 */

// TODO this is only needed in slip_ref_triangular_solve and doesn't need to
// be its own function.

#include "SLIP_LU_internal.h"

SLIP_info slip_reset_mpz_array
(
    mpz_t *x,      // mpz array to be reset
    int32_t n,     // size of x
    int32_t top,   // beginning of the nonzero pattern
    int32_t *xi    // nonzero pattern
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info;
    if (!x || n <= 0 || top < 0 || !xi) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------
    // clear x
    //--------------------------------------------------------------------------

    // Access the nonzero pattern located in xi[top..n-1]
    // and set nonzero x[i] = 0
    for (int32_t i = top; i < n; i++)
    {
        SLIP_CHECK (SLIP_mpz_set_ui (x [xi [i]], 0)) ;
    }

    return SLIP_OK;
}
