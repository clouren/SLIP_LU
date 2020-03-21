//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_output_soln: Output x to matlab
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function outputs the solution as a mxArray. */

#include "SLIP_LU_mex.h"

mxArray* slip_mex_output_soln
(
    double** x,        // The matrix to be output
    int64_t m,         // size of x
    int64_t n          // the size of x
)
{
    // Create a m * n array
    mxArray* Xmatlab = mxCreateDoubleMatrix ((mwSize) m, (mwSize) n, mxREAL);

    // Get the numeric values
    double* x2 = mxGetPr(Xmatlab);
    int64_t count = 0;

    // Populate the nonzeros in output matrix
    for (int64_t j = 0; j < n; j++)
    {
        for (int64_t i = 0; i < m; i++)
        {
            x2[count++] = x[i][j];
        }
    }
    return Xmatlab;
}

