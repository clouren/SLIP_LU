//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_check_for_inf: Check A and b for inf
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"


/* Purpose: This function checks if the input matrix or RHS has numbers too
 * large for double*/
void slip_mex_check_for_inf
(
    double* x, // The array of numeric values 
    mwSize n   // size of array
)
{
    for (mwSize k = 0; k < n; k++)
    {
        if (mxIsInf(x[k]))
        {
            mexErrMsgTxt("Numbers are too large for double. Please try the C "
                "code with mpfr, mpq, or mpz");
        }
    }
}
