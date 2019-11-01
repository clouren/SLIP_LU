#include "SLIP_LU_mex.h"

/* This is one of the 4 cpp files which defines the SLIP LU matlab interface
 * This one defines [L U P Q x] = SLIP_LU(A,b) where L and U are doolittle and
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
    mp_set_memory_functions(mxMalloc, mxGMPRealloc, mxGMPFree);
    //SLIP_initialize();
    SLIP_info status;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    SLIP_check_input(pargin, nargin);
    if (nargout != 5 || nargin != 3)
    {
        mexErrMsgTxt ("Usage: [L U P Q x] = SLIP_LU(A,b,option)") ;
    }

    //--------------------------------------------------------------------------
    // Allocate and initialize data structures
    //--------------------------------------------------------------------------
    SLIP_sparse *A = NULL, *L = NULL, *U = NULL;
    A = SLIP_create_sparse();
    L = SLIP_create_sparse();
    U = SLIP_create_sparse();
    SLIP_dense *b = SLIP_create_dense();
    //Set defaults for options
    SLIP_options* option = SLIP_create_default_options();
    if (!A || !L || !U || !b || !option)
    {
        SLIP_mex_error (SLIP_OUT_OF_MEMORY);
    }
    SLIP_LU_analysis* S = NULL;

    //--------------------------------------------------------------------------
    // Process input and declare variables
    //--------------------------------------------------------------------------
    // Read in options
    SLIP_get_matlab_options(option, pargin[2]);

    // Read in A and b
    SLIP_mex_get_A_and_b(A, b, pargin, nargin);

    // Create arrays based on the size of input matrix
    S = SLIP_create_LU_analysis((A->n)+1);
    double** soln = SLIP_create_double_mat(b->m, b->n);
    int32_t* pinv = (int32_t*) SLIP_malloc(A->n* sizeof(int32_t));
    mpz_t* rhos = SLIP_create_mpz_array(A->n);
    mpq_t** soln_mpq = SLIP_create_mpq_mat(b->m, b->n);
    if (!S || !soln || !pinv || !rhos || !soln_mpq)
    {
        SLIP_mex_error (SLIP_OUT_OF_MEMORY);
    }

    //--------------------------------------------------------------------------
    // Symbolic analysis, factorization, and solve
    //--------------------------------------------------------------------------
    // Symbolic Analysis
    SLIP_MEX_OK (SLIP_LU_analyze(S, A, option));
    // SLIP LU Factorization
    SLIP_MEX_OK (SLIP_LU_factorize(L, U, A, S, rhos, pinv, option));
    
    
    // FB Sub
    
    SLIP_MEX_OK (SLIP_LU_solve(soln_mpq, b, rhos, L, U, pinv));
    
    SLIP_MEX_OK (SLIP_permute_x(soln_mpq, b->m, b->n, S));
    
    SLIP_MEX_OK (SLIP_scale_x(soln_mpq, A, b));
    SLIP_MEX_OK (SLIP_get_double_soln(soln, soln_mpq, b->m, b->n));
    
    //--------------------------------------------------------------------------
    // Output
    //--------------------------------------------------------------------------
    pargout[0] = SLIP_mex_output_L(L, rhos);                     // out[0] = L
    pargout[1] = SLIP_mex_output_U(U, rhos, A->scale);           // out[1] = U
    pargout[2] = SLIP_mex_output_p(pinv, A->n);			 // out[2] = P
    pargout[3] = SLIP_mex_output_col_permut(S->q, A->n);         // out[3] = Q
    pargout[4] = SLIP_mex_output_soln(soln, b->m, b->n);         // out[4] = soln

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    SLIP_delete_mpq_mat(&soln_mpq, b->m, b->n);
    SLIP_delete_mpz_array(&rhos, A->n);
    SLIP_FREE(pinv);
    SLIP_delete_double_mat(&soln, b->m, b->n);
    SLIP_delete_LU_analysis(&S);
    SLIP_FREE(option);
    SLIP_delete_dense(&b);
    SLIP_delete_sparse(&U);
    SLIP_delete_sparse(&L);
    SLIP_delete_sparse(&A);
    //SLIP_finalize();
}
