//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_output_col_permut: Output column perm. to MATLAB
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function outputs an int64_t array as a mxArray.*/

#include "SLIP_LU_mex.h"

mxArray* slip_mex_output_col_permut
(
    int64_t* x,        // int64_t array to be output
    int64_t n          // size of x
)
{
    // Create a n*1 matlab array
    mxArray* Xmatlab = mxCreateDoubleMatrix ((mwSize) n, 1, mxREAL);

    // Numeric values of Xmatlab
    double* x2 = mxGetPr(Xmatlab);

    // Get cast x as double
    for (int64_t k = 0; k < n; k++)
    {
        x2[k] = (double) x[k];
    }
    return Xmatlab;
}

