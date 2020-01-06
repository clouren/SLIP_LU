//------------------------------------------------------------------------------
// SLIP_LU/SLIP_mex_soln3: Interface to SLIP LU via matlab
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"
/* Purpose: One of the 3 c files defining the SLIP LU matlab intreface
 * This one defines [L U P Q] = SLIP_LU(A) where L and U are Doolittle and
 * L*U = P*A*Q
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
    SLIP_initialize_expert(mxMalloc, SLIP_gmp_mex_realloc, SLIP_gmp_mex_free);

    //mp_set_memory_functions(mxMalloc, SLIP_gmp_mex_realloc, SLIP_gmp_mex_free);
    //SLIP_initialize();
    SLIP_info status;

    //--------------------------------------------------------------------------
    // Input Check
    //--------------------------------------------------------------------------
    slip_check_input(pargin, nargin);
    if (nargout != 4 || nargin != 2)
    {
        mexErrMsgTxt ("Usage: [L U P Q] = SLIP_LU(A,option)") ;
    }

    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------
    SLIP_sparse *A = NULL, *L = NULL, *U = NULL;    // Input matrix
    A = SLIP_create_sparse();
    L = SLIP_create_sparse();
    U = SLIP_create_sparse();
    //Set defaults for options
    SLIP_options* option = SLIP_create_default_options();
    if (!A || !L || !U || !option)
    {
        slip_mex_error (SLIP_OUT_OF_MEMORY);
    }
    SLIP_LU_analysis* S = NULL;

    //--------------------------------------------------------------------------
    // Declare variables and process input
    //--------------------------------------------------------------------------
    // Read in options
    slip_get_matlab_options(option, pargin[1]);

    // Read in A and b
    slip_mex_get_A_and_b(A, NULL, pargin, nargin);

    S = SLIP_create_LU_analysis((A->n)+1);
    int32_t* pinv = (int32_t*) SLIP_malloc((A->n)* sizeof(int32_t));
    mpz_t* rhos = SLIP_create_mpz_array(A->n);
    if (!S || !pinv || !rhos)
    {
        slip_mex_error (SLIP_OUT_OF_MEMORY);
    }

    //--------------------------------------------------------------------------
    // Symbolic analysis and factorization
    //--------------------------------------------------------------------------
    //Symbolic Analysis
    SLIP_MEX_OK (SLIP_LU_analyze(S, A, option));
    // SLIP LU Factorization
    SLIP_MEX_OK (SLIP_LU_factorize(L, U, A, S, rhos, pinv, option));

    //--------------------------------------------------------------------------
    // Output
    //--------------------------------------------------------------------------
    
    pargout[0] = slip_mex_output_L(L, rhos);                       // out[0] = L
    pargout[1] = slip_mex_output_U(U, rhos, A->scale);             // out[1] = U
    pargout[2] = slip_mex_output_p(pinv, A->n);                    // out[2] = P
    pargout[3] = slip_mex_output_col_permut(S->q, A->n);           // out[3] = Q

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
