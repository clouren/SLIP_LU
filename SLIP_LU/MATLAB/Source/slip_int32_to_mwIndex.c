//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_int64_to_mwindex: Convert int64_t* array to mwIndex*
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts an int64_t to an mwIndex */

#include "SLIP_LU_mex.h"

void slip_int64_to_mwIndex
(
    mwIndex* y,    // the mwIndex vector
    int64_t* x,    // the int64_t vector
    int64_t n      // the size of x and y
)
{
    for (int64_t i = 0; i < n; i++)
    {
        y[i] = (mwIndex) x[i];
    }
}

