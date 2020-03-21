//------------------------------------------------------------------------------
// SLIP_LU/slip_reset_int64_array: clear all of an int64_t workspace array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function initializes an int64_t vector of size n and sets the
 * value equal to -1. This function is used for the history and pivot vectors.
 */

#include "SLIP_LU_internal.h"

SLIP_info slip_reset_int64_array
(
    int64_t *h,    // int64_t vector to be reset
    int64_t n      // size of the int64_t vector
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!h || n <= 0) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------
    // Update h
    //--------------------------------------------------------------------------

    for (int64_t i = 0; i < n; i++)
    {
        h[i] = -1;
    }
    return SLIP_OK;
}
