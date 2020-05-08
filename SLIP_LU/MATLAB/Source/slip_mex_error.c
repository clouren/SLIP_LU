//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_error: Return error messages to matlab
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function prints error messages for MATLAB for debugging*/

#include "SLIP_LU_mex.h"

void slip_mex_error
(
    SLIP_info status
)
{
    
    switch (status)
    {
        case SLIP_OK :                   // all is well
            return ;

        case SLIP_OUT_OF_MEMORY :        // out of memory
            mexErrMsgTxt ("out of memory") ;

        case SLIP_SINGULAR :             // the input matrix A is singular
            mexErrMsgTxt ("input matrix is singular") ;

        case SLIP_INCORRECT_INPUT :      // one or more input arguments are incorrect
            mexErrMsgTxt ("invalid inputs") ;

        case SLIP_INCORRECT :            // The solution is incorrect
            mexErrMsgTxt ("result invalid") ;

        case SLIP_PANIC :                // SLIP_LU used without proper initialization
        default : 
            mexErrMsgTxt ("unknown error") ;
    }
}

