//------------------------------------------------------------------------------
// SLIP_LU/slip_reset_int32_array2: clear the used parts of an int32_t workspace
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function resets an int32_t vector of size n and sets each term
 * equal to -1 with nonzero pattern given. This is more efficient than
 * resetting each term individually.
 */

// TODO this is only needed in slip_ref_triangular_solve and doesn't need to
// be its own function.

#include "SLIP_LU_internal.h"

SLIP_info slip_reset_int32_array2
(
    int32_t *h,    // int32_t vector to be reset
    int32_t n,     // size of h
    int32_t top,   // beginning of nonzero pattern
    int32_t *xi    // nonzero pattern
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!h || !xi || n <= 0 || top < 0) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------
    // clear the h array
    //--------------------------------------------------------------------------

    // Access the nonzero pattern located in xi [top..n-1]
    // and set each h[i] = -1, for each i in the list xi [top..n-1].

    for (int32_t i = top; i < n; i++)
    {
        h[xi[i]] = -1;
    }
    return SLIP_OK;
}
