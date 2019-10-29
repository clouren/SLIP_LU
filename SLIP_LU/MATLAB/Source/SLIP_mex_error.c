#include "SLIP_LU_mex.h"

/* Purpose: This function prints error messages for MATLAB */
void SLIP_mex_error
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
