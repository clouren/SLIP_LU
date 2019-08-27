#include "SLIP_LU_mex.h"

/* Purpose: This function outputs the p matrix from pinv as a mxArray. */
mxArray* SLIP_mex_output_p
(
    int32_t* pinv,     // pinv
    int32_t n          // size of pinv 
)
{
    // Create a n*1 array
    mxArray* Pmatlab = mxCreateDoubleMatrix ((mwSize) n, 1, mxREAL);

    // Numeric values of Pmatlab
    double* x = mxGetPr(Pmatlab);

    // Set Pmatlab[pinv[k]] = k
    for (int32_t k = 0; k < n; k++)
    {
        x[pinv[k]] = (double) k;
    }
    return Pmatlab;
}
