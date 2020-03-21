//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_output_p: Output the row perm to MATLAB
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function outputs the p matrix from pinv as a mxArray. */

#include "SLIP_LU_mex.h"

mxArray* slip_mex_output_p
(
    int64_t* pinv,     // pinv
    int64_t n          // size of pinv
)
{
    // Create a n*1 array
    mxArray* Pmatlab = mxCreateDoubleMatrix ((mwSize) n, 1, mxREAL);

    // Numeric values of Pmatlab
    double* x = mxGetPr(Pmatlab);

    // Set Pmatlab[pinv[k]] = k
    for (int64_t k = 0; k < n; k++)
    {
        x[pinv[k]] = (double) k;
    }
    return Pmatlab;
}

