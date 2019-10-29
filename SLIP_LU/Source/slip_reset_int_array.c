//------------------------------------------------------------------------------
// SLIP_LU/slip_reset_int_array: clear all of an int workspace array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function initializes an int vector of size n and sets the value
 * equal to -1. This function is used for the history and pivot vectors. 
 */
SLIP_info slip_reset_int_array
(
    int32_t *h,    // int32_t vector to be reset
    int32_t n      // size of the int32_t vector
)    
{
    // Check input
    if (!h || n <= 0) {return SLIP_INCORRECT_INPUT;}
    // Update h
    for (int32_t i = 0; i < n; i++)
    {
        h[i] = -1;
    }
    return SLIP_OK;
}
