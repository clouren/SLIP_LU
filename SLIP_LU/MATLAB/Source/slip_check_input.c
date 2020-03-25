//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_check_input: Check the input obtained from MATLAB user
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"

/* Purpose: This function checks the input of the MATLAB array. It ensures that
   the input matrix and right hand side vectors are in the correct format. */

void slip_check_input
(
    const mxArray * input [],   // The MATLAB inputs
    int nargin                  // # of input arguments
)
{

    if (nargin != 2 && nargin != 3)
    {
        mexErrMsgTxt ("incorrect number of inputs") ;
    }

    //--------------------------------------------------------------------------
    // Check matrix A
    //--------------------------------------------------------------------------

    if (mxIsComplex (input[0]))        // Is the matrix real valued?
    {
        mexErrMsgTxt ("Matrix must be real; try backslash instead") ;
    }
    if (!mxIsSparse (input[0]))     // Is the matrix sparse?
    {
        mexErrMsgTxt ("Matrix must be sparse. Try again with A = sparse(A)") ;
    }

    //--------------------------------------------------------------------------
    // Check option
    //--------------------------------------------------------------------------

    if (!mxIsStruct (input[nargin-1]))        // Is third argument the struct?
    {
        mexErrMsgTxt ("Third argument must be the options struct") ;
    }

    //--------------------------------------------------------------------------
    // Check option
    //--------------------------------------------------------------------------

    if (nargin == 3)
    {
        if (mxIsSparse (input[1]))         // Is b sparse?
        {
            mexErrMsgTxt ("Right hand side vector must be dense. "
                "Try again with b = full(b)") ;
        }
    }
}

