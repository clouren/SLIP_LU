//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mwIndex_to_int32: Convert mwIndex* to int32_t*
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"

/* Purpose: This function converts an mwIndex to an int32_t*/
void slip_mwIndex_to_int32
(
    int32_t* y,    // the int32_t vector
    mwIndex* x,    // the mwIndex vector
    mwSize n       // the size of x and y
)
{
    for (mwSize i = 0; i < n; i++)
    {  
        y[i] = (int32_t) x[i];
    }
}
