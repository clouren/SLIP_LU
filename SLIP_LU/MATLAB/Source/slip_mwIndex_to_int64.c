//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mwIndex_to_int64: Convert mwIndex* to int64_t*
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts an mwIndex to an int64_t, used for A&b*/

#include "SLIP_LU_mex.h"

void slip_mwIndex_to_int64
(
    int64_t* y,    // the int64_t vector
    mwIndex* x,    // the mwIndex vector
    mwSize n       // the size of x and y
)
{
    for (mwSize i = 0; i < n; i++)
    {
        y[i] = (int64_t) x[i];
    }
}

