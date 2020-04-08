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

    SLIP_initialize_expert (mxMalloc, mxCalloc, mxRealloc, mxFree) ;
    SuiteSparse_config.printf_func = mexPrintf ;
    SLIP_info status ;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    slip_check_input(pargin, nargin);
    if (nargout > 1 || nargout <= 0 || nargin != 3)
    {
        mexErrMsgTxt("Usage: x = SLIP_LU(A,b,option)");
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
    // Read in options
    slip_get_matlab_options(option, pargin[2]);

    // Read in A and b
    slip_mex_get_A_and_b(&A, &b, pargin, nargin, option);

    //--------------------------------------------------------------------------
    // Solve
    //--------------------------------------------------------------------------

    SLIP_MEX_OK(SLIP_backslash( &x, SLIP_FP64, A, b, option));
    
    //--------------------------------------------------------------------------
    // Set outputs, free memory
    //--------------------------------------------------------------------------

    mxArray* Xmatlab = mxCreateDoubleMatrix ((mwSize) x->m, (mwSize) x->n,
        mxREAL);
    double* x2 = mxGetPr(Xmatlab);

    for (int k = 0; k < x->n*x->m; k++)
    {
        x2[k] = x->x.fp64[k];
    }

    //TODO Can this work?  Tim: no, but mxSetData* can do it.  FIXME
    //x2 = x->x.fp64;

    pargout[0] =  Xmatlab;
    SLIP_matrix_free(&b, option);
    SLIP_matrix_free(&A, option);
    SLIP_FREE(option);
    SLIP_finalize();
}

