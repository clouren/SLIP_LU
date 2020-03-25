//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/SLIP_gmp_mex_free: A gmp free function for Matlab mex files
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: A GMP wrapper for mxFree
 * MATLAB's mxFree instead of free */

#include "SLIP_LU_mex.h"

// A GMP realloc function
void SLIP_gmp_mex_free
(
    void* x,    // void* to be freed
    size_t a    // Size
)
{
    if (x)
    {
        mxFree(x);
    }
}

