#include "SLIP_LU_mex.h"

/* Purpose: This function outputs an int32_t array as a mxArray.*/
mxArray* SLIP_mex_output_col_permut
(
    int32_t* x,        // int32_t array to be output
    int32_t n          // size of x
)
{
    // Create a n*1 matlab array
    mxArray* Xmatlab = mxCreateDoubleMatrix ((mwSize) n, 1, mxREAL);
    
    // Numeric values of Xmatlab
    double* x2 = mxGetPr(Xmatlab);
    
    // Get cast x as double
    for (int32_t k = 0; k < n; k++)
    {
        x2[k] = (double) x[k];
    }
    return Xmatlab;
}
