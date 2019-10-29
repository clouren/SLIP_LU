#include "SLIP_LU_mex.h"
/* This is one of the cpp files which defines SLIP LU Matlab interface
 * This one defines d = SLIP_det(A)
 */


void mexFunction
(
    int32_t nargout,
    mxArray *pargout [ ],
    int32_t nargin,
    const mxArray *pargin [ ]
)
{
    //--------------------------------------------------------------------------
    // Initialize SLIP LU library environment
    //--------------------------------------------------------------------------
    SLIP_initialize();
    SLIP_info status;

    //--------------------------------------------------------------------------
    // Input checks
    //--------------------------------------------------------------------------
    SLIP_check_input(pargin, nargin);
    if (nargout != 1 || nargin != 2)
    {
        mexErrMsgTxt ("Usage: d = SLIP_det(A)") ;
    }

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------
    SLIP_sparse *A = NULL, *L = NULL, *U = NULL;                // Input matrix
    A = SLIP_create_sparse();
    L = SLIP_create_sparse();
    U = SLIP_create_sparse();
    //Set defaults for options
    SLIP_options* option = SLIP_create_default_options();
    if (!A || !L || !U || !option)
    {
        SLIP_mex_error (SLIP_OUT_OF_MEMORY);
    }
    SLIP_LU_analysis* S = NULL;

    //--------------------------------------------------------------------------
    // Declare Variables & process input
    //--------------------------------------------------------------------------
    SLIP_get_matlab_options(option, pargin[1]); // Read in options

    // Read in A and b
    SLIP_mex_get_A_and_b(A, NULL, pargin, nargin);

    S = SLIP_create_LU_analysis((A->n)+1);
    int32_t* pinv = (int32_t*) SLIP_malloc(A->n* sizeof(int32_t));
    mpz_t* rhos = SLIP_create_mpz_array(A->n);
    if (!S || !pinv || !rhos)
    {
        SLIP_mex_error (SLIP_OUT_OF_MEMORY);
    }

    //--------------------------------------------------------------------------
    // Symbolic Analysis, factorization, and solve
    //--------------------------------------------------------------------------
    SLIP_MEX_OK (SLIP_LU_analyze(S, A, NULL, option));//Symbolic Analysis
    SLIP_MEX_OK (SLIP_LU_factorize(L, U, A, S, rhos, pinv, option));

    //--------------------------------------------------------------------------
    // Output determinant
    //--------------------------------------------------------------------------
    pargout[0] = SLIP_mex_output_determinant(rhos[(A->n)-1], A);

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    SLIP_delete_mpz_array(&rhos, A->n);
    SLIP_FREE(pinv);
    SLIP_delete_LU_analysis(&S);
    SLIP_FREE(option);
    SLIP_delete_sparse(&U);
    SLIP_delete_sparse(&L);
    SLIP_delete_sparse(&A);
    SLIP_finalize();
}
