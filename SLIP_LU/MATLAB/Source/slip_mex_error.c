//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_error: Return error messages to matlab
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"

/* Purpose: This function prints error messages for MATLAB */
void slip_mex_error
(
    SLIP_info status
)
{
    if (status == SLIP_OUT_OF_MEMORY)
    {
        mexErrMsgTxt("Error, Out of memory");
    }
    else if (status == SLIP_SINGULAR)
    {
        mexErrMsgTxt("Error, Input matrix is singular");
    }
    else if (status == SLIP_INCORRECT_INPUT)
    {
        mexErrMsgTxt("Error, Input is incorrect");
    }
}
