//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/SLIP_mex_soln: Use SLIP LU within MATLAB
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: The .c file defining the SLIP LU matlab interfacee
 * This function defines: x = SLIP_LU(A, b, option)
 */

#include "SLIP_LU_mex.h"

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    //--------------------------------------------------------------------------
    // Initialize SLIP LU library environment
    //--------------------------------------------------------------------------

    SLIP_info status ;
    SLIP_MEX_OK (SLIP_initialize_expert
        (mxMalloc, mxCalloc, mxRealloc, mxFree)) ;
    SuiteSparse_config.printf_func = mexPrintf ;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    if (nargout > 1 || nargin < 2 || nargin > 3)
    {
        mexErrMsgTxt("Usage: x = SLIP_mex_soln (A,b,option)");
    }

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (mxIsComplex (pargin [0]) || mxIsComplex (pargin [1]))
    {
        mexErrMsgTxt ("Inputs must be real") ;
    }
    if (!mxIsSparse (pargin [0]))     // Is the matrix sparse?
    {
        mexErrMsgTxt ("First input must be sparse") ;
    }
    if (mxIsSparse (pargin [1]))         // Is b sparse?
    {
        mexErrMsgTxt ("Second input must be full") ;
    }

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------

    SLIP_matrix *A = NULL;
    SLIP_matrix *b = NULL;
    SLIP_matrix *x = NULL;
    SLIP_options *option = SLIP_create_default_options();
    if (option == NULL)
    {
        slip_mex_error (SLIP_OUT_OF_MEMORY);
    }

    //--------------------------------------------------------------------------
    // Declare variables and process input
    //--------------------------------------------------------------------------

    // get options
    if (nargin > 2) slip_get_matlab_options (option, pargin [2]) ;

    // convert MATLAB inputs A and b into SLIP_matrix objects
    slip_mex_get_A_and_b (&A, &b, pargin, nargin, option) ;

    //--------------------------------------------------------------------------
    // x = A\b via SLIP_LU, returning result as double
    //--------------------------------------------------------------------------

    SLIP_MEX_OK (SLIP_backslash (&x, SLIP_FP64, A, b, option)) ;
    
    //--------------------------------------------------------------------------
    // Set outputs, free memory
    //--------------------------------------------------------------------------

    // create an empty 0-by-0 MATLAB matrix and free its contents
    pargout [0] = mxCreateDoubleMatrix (0, 0, mxREAL) ;
    mxFree (mxGetDoubles (pargout [0])) ;

    // transplant x into the new MATLAB matrix and set its size
    mxSetDoubles (pargout [0], x->x.fp64) ;
    mxSetM (pargout [0], x->m) ;
    mxSetN (pargout [0], x->n) ;
    x->x.fp64 = NULL ;  // set to NULL so it is not freed by SLIP_matrix_free

    SLIP_matrix_free (&x, option) ;
    SLIP_matrix_free (&b, option) ;
    SLIP_matrix_free (&A, option) ;
    SLIP_FREE (option) ;
    SLIP_finalize ( ) ;
}

