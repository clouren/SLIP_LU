//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_get_A_and_b: Obtain user's A and b matrices
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function reads in the A matrix and right hand side vectors. */

#include "SLIP_LU_mex.h"

void slip_mex_get_A_and_b
(
    SLIP_matrix **A_handle,  // Internal SLIP Mat stored in CSC
    SLIP_matrix **b_handle,   // mpz matrix used internally
    const mxArray* pargin[], // The input A matrix and options
    int nargin,               // Number of input to the mexFunction
    SLIP_options* option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!A_handle || !pargin)
    {
        slip_mex_error (SLIP_INCORRECT_INPUT);
    }
    (*A_handle) = NULL ;

    //--------------------------------------------------------------------------
    // Declare variables
    //--------------------------------------------------------------------------

    SLIP_info status;
    mxArray* tmp;
    mwSize nA, mA, nb, mb, Anz, k, j;
    mwIndex *Ap, *Ai;
    double *Ax, *bx;

    //--------------------------------------------------------------------------
    // Read in A
    //--------------------------------------------------------------------------

    // Read in Ap, Ai, Ax
    Ap = mxGetJc (pargin[0]) ;
    Ai = mxGetIr (pargin[0]) ;
    Ax = mxGetDoubles (pargin[0]) ;
    if (!Ai || !Ap || !Ax)
    {
        slip_mex_error (SLIP_INCORRECT_INPUT) ;
    }

    // Get info about A
    nA = mxGetN (pargin[0]) ;
    mA = mxGetM (pargin[0]) ;
    Anz = (mwSize) Ap[nA];
    if (nA != mA)
    {
        mexErrMsgTxt ("A must be square") ;
    }

    int64_t *Ai_int64 = (int64_t *) SLIP_malloc (Anz * sizeof (int64_t)) ;
    int64_t *Ap_int64 = (int64_t *) SLIP_malloc ((nA+1) * sizeof (int64_t)) ;
    if (!Ai_int64 || !Ap_int64)
    {
        slip_mex_error (SLIP_OUT_OF_MEMORY) ;
    }

    // convert the integer pattern
    slip_mwIndex_to_int64 (Ai_int64, Ai, Anz) ;
    slip_mwIndex_to_int64 (Ap_int64, Ap, nA+1) ;

    // check the values of A
    bool A_has_int64_values = slip_mex_check_for_inf (Ax, Anz) ;

    SLIP_matrix* A = NULL;
    SLIP_matrix* A_matlab = NULL;
    
    if (A_has_int64_values)
    {
        // All entries in A can be typecast to int64_t without change in value.
        int64_t *Ax_int64 = (int64_t*) SLIP_malloc (Anz* sizeof (int64_t)) ;
        if (!Ax_int64)
        {
            slip_mex_error (SLIP_OUT_OF_MEMORY) ;
        }
        for (k = 0; k < Anz; k++)
        {
            Ax_int64[k] = (int64_t) Ax[k];
        }
        
        
        // Create A_matlab (which is shallow)
        SLIP_matrix_allocate(&A_matlab, SLIP_CSC, SLIP_INT64, (int64_t) mA, (int64_t) nA,
                             (int64_t) Anz, true, true, option);
        
        A_matlab->p = Ap_int64;
        A_matlab->i = Ai_int64;
        A_matlab->x.int64 = Ax_int64;
    }
    else
    {
        // Entries in A cannot be typecast to int64_t without changing them.
        
        // Create A_matlab
        SLIP_matrix_allocate(&A_matlab, SLIP_CSC, SLIP_FP64, (int64_t) mA, (int64_t) nA,
                             (int64_t) Anz, true, true, option);
        
        A_matlab->p = Ap_int64;
        A_matlab->i = Ai_int64;
        A_matlab->x.fp64 = Ax;
    }

    // Create A
    status = SLIP_matrix_copy(&A, SLIP_CSC, SLIP_MPZ, A_matlab, option);

    if (status != SLIP_OK)
    {
        mexErrMsgTxt ("A invalid") ;
    }
    SLIP_FREE (Ai_int64) ;
    SLIP_FREE (Ap_int64) ;

    //--------------------------------------------------------------------------
    // Read in b
    //--------------------------------------------------------------------------

    SLIP_matrix* b = NULL;
    SLIP_matrix* b_matlab = NULL;
    
    if (nargin == 3)
    {
        bx = mxGetDoubles (pargin[1]) ;
        if (!bx)
        {
            slip_mex_error (SLIP_INCORRECT_INPUT) ;
        }

        // Get info about RHS vector (s)
        nb = mxGetN (pargin[1]) ;
        mb = mxGetM (pargin[1]) ;
        if (mb != mA)
        {
            mexErrMsgTxt ("dimension mismatch") ;
        }

        int64_t count = 0;

        // check the values of b
        bool b_has_int64_values = slip_mex_check_for_inf (bx, nb*mb) ;

        if (b_has_int64_values)
        {
            
            // Create b_matlab (which is shallow)
            SLIP_matrix_allocate(&b_matlab, SLIP_DENSE, SLIP_INT64, (int64_t) mb,
                                 (int64_t) nb, (int64_t) mb*nb, true, true, option);
        
            b_matlab->x.int64 = SLIP_calloc( (int64_t) nb*mb, sizeof(int64_t));
            for (int64_t j = 0; j < mb*nb; j++)
            {
                b_matlab->x.int64[j] = (int64_t) bx[j];
            }
        
        }
        else
        {
            
            // Create b_matlab (which is shallow)
            SLIP_matrix_allocate(&b_matlab, SLIP_DENSE, SLIP_FP64, (int64_t) mb,
                                 (int64_t) nb, (int64_t) mb*nb, true, true, option);
        
            b_matlab->x.fp64 = bx;
        }

        // Create b
        status = SLIP_matrix_copy(&b, SLIP_DENSE, SLIP_MPZ, b_matlab, option);

        if (status != SLIP_OK)
        {
            mexErrMsgTxt ("b invalid") ;
        }
    }

    (*A_handle) = A;
    (*b_handle) = b;
}

