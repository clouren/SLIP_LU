#include "SLIP_LU_mex.h"


/* Purpose: This function checks if the input matrix or RHS has numbers too
 * large for double*/
void SLIP_mex_check_for_inf
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
