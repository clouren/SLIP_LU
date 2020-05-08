//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_get_matlab_options: Set factorization options for SLIP LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"

// Purpose: This function reads in the necessary information from the options
// struct for MATLAB.

#define MATCH(s,t) (strcmp (s,t) == 0)

void slip_get_matlab_options
(
    SLIP_options* option,   // Control parameters (must not be NULL)
    const mxArray* input    // options struct, may be NULL
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    mxArray *field ;
    #define LEN 256
    char string [LEN+1] ;

    // true if input options struct is present 
    bool present = (input != NULL) && !mxIsEmpty (input) && mxIsStruct (input) ;

    //--------------------------------------------------------------------------
    // Get the column ordering
    //--------------------------------------------------------------------------

    option->order = SLIP_COLAMD ;     // default: COLAMD ordering
    field = present ? mxGetField (input, 0, "order") : NULL ;
    if (field != NULL)
    {
        if (!mxIsChar (field)) mexErrMsgTxt ("option.order must be a string") ;
        mxGetString (field, string, LEN) ;
        if (MATCH (string, "none"))
        {
            option->order = SLIP_NO_ORDERING ;  // None: A is factorized as-is
        }
        else if (MATCH (string, "colamd"))
        {
            option->order = SLIP_COLAMD ;       // COLAMD: Default
        }
        else if (MATCH (string, "amd"))
        {
            option->order = SLIP_AMD ;          // AMD
        }
        else
        {
            mexErrMsgTxt ("unknown option.order") ;
        }
    }

    //--------------------------------------------------------------------------
    // Get the row pivoting scheme
    //--------------------------------------------------------------------------

    option->pivot = SLIP_TOL_SMALLEST ;     // default: diag pivoting with tol
    field = present ? mxGetField (input, 0, "pivot") : NULL ;
    if (field != NULL)
    {
        if (!mxIsChar (field)) mexErrMsgTxt ("option.pivot must be a string") ;
        mxGetString (field, string, LEN) ;
        if (MATCH (string, "smallest"))
        {
            option->pivot = SLIP_SMALLEST ;         // Smallest pivot
        }
        else if (MATCH (string, "diagonal"))
        {
            option->pivot = SLIP_DIAGONAL ;         // Diagonal pivoting
        }
        else if (MATCH (string, "first"))
        {
            option->pivot = SLIP_FIRST_NONZERO ;    // First nonzero in column
        }
        else if (MATCH (string, "tol smallest"))
        {
            option->pivot = SLIP_TOL_SMALLEST ;     // diag pivoting with tol
                                                    // for smallest pivot
        }
        else if (MATCH (string, "tol largest"))
        {
            option->pivot = SLIP_TOL_LARGEST ;      // diag pivoting with tol
                                                    // for largest pivot.
        }
        else if (MATCH (string, "largest"))
        {
            option->pivot = SLIP_LARGEST ;          // Largest pivot
        }
        else
        {
            mexErrMsgTxt ("unknown option.pivot") ;
        }
    }

    //--------------------------------------------------------------------------
    // tolerance for row partial pivoting
    //--------------------------------------------------------------------------

    option->tol = 0.1 ;     // default tolerance is 0.1
    if (option->pivot == SLIP_TOL_SMALLEST || option->pivot == SLIP_TOL_LARGEST)
    {
        field = present ? mxGetField (input, 0, "tol") : NULL ;
        if (field != NULL)
        {
            option->tol = mxGetScalar (field) ;
            if (option->tol > 1 || option->tol <= 0)
            {
                mexErrMsgTxt ("invalid option.tol, must be > 0 and <= 1") ;
            }
        }
    }
}

