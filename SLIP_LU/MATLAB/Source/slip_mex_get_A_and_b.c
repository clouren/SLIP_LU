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
    SLIP_sparse **A_handle,  // Internal SLIP Mat stored in CSC
    SLIP_dense **b_handle,   // mpz matrix used internally
    const mxArray* pargin[], // The input A matrix and options
    int32_t nargin           // Number of input to the mexFunction
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
    // TODO allow option to be NULL in all C functions
    SLIP_options *option = SLIP_create_default_options ();

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

    int32_t *Ai_int32 = (int32_t *) SLIP_malloc (Anz * sizeof (int32_t)) ;
    int32_t *Ap_int32 = (int32_t *) SLIP_malloc ((nA+1) * sizeof (int32_t)) ;
    if (!Ai_int32 || !Ap_int32)
    {
        slip_mex_error (SLIP_OUT_OF_MEMORY) ;
    }

    // convert the integer pattern
    slip_mwIndex_to_int32 (Ai_int32, Ai, Anz) ;
    slip_mwIndex_to_int32 (Ap_int32, Ap, nA+1) ;

    // check the values of A
    bool A_has_int32_values = slip_mex_check_for_inf (Ax, Anz) ;

    if (A_has_int32_values)
    {
        // All entries in A can be typecast to int32_t without change in value.
        int32_t *Ax_int32 = (int32_t*) SLIP_malloc (Anz* sizeof (int32_t)) ;
        if (!Ax_int32)
        {
            slip_mex_error (SLIP_OUT_OF_MEMORY) ;
        }
        for (k = 0; k < Anz; k++)
        {
            Ax_int32[k] = (int32_t) Ax[k];
        }
        // Create A with no scaling
        status = SLIP_build_sparse_csc_int32 (A_handle, Ap_int32, Ai_int32,
            Ax_int32, (int32_t) nA, (int32_t) Anz) ;
        SLIP_FREE (Ax_int32) ;
    }
    else
    {
        // Entries in A cannot be typecast to int32_t without changing them.
        // Create A with scaling
        status = SLIP_build_sparse_csc_double (A_handle, Ap_int32, Ai_int32,
            Ax, (int32_t) nA, (int32_t) Anz, option) ;
    }

    if (status != SLIP_OK)
    {
        mexErrMsgTxt ("A invalid") ;
    }
    SLIP_FREE (Ai_int32) ;
    SLIP_FREE (Ap_int32) ;

    //--------------------------------------------------------------------------
    // Read in b
    //--------------------------------------------------------------------------

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

        int32_t count = 0;

        // check the values of b
        bool b_has_int32_values = slip_mex_check_for_inf (bx, nb*mb) ;

        if (b_has_int32_values)
        {
            // populate bx to a int32_t mat
            int32_t** bx_int32 = SLIP_create_int32_mat ((int32_t) mb,
                (int32_t) nb) ;
            if (!bx_int32)
            {
                slip_mex_error (SLIP_OUT_OF_MEMORY) ;
            }
            for (int32_t j = 0; j < nb; j++)
            {
                for (int32_t i = 0; i < mb; i++)
                {
                    bx_int32[i][j] = (int32_t) bx[count];
                    count++;
                }
            }

            // Create b
            status = SLIP_build_dense_int32(b_handle, bx_int32, (int32_t) mb,
                (int32_t) nb) ;
            SLIP_delete_int32_mat (&bx_int32, (int32_t) mb, (int32_t) nb) ;
        }
        else
        {
            // populate bx to a double mat
            double** bx_doub;
            bx_doub = SLIP_create_double_mat ((int32_t) mb, (int32_t) nb) ;
            if (!bx_doub)
            {
                slip_mex_error (SLIP_OUT_OF_MEMORY) ;
            }

            count = 0;
            for (int32_t j = 0; j < nb; j++)
            {
                for (int32_t i = 0; i < mb; i++)
                {
                    bx_doub[i][j] = bx[count];
                    count++;
                }
            }

            // Create b
            status = SLIP_build_dense_double (b_handle, bx_doub, (int32_t) mb,
                (int32_t) nb, option) ;
            SLIP_delete_double_mat (&bx_doub, (int32_t) mb, (int32_t) nb) ;
        }

        if (status != SLIP_OK)
        {
            mexErrMsgTxt ("b invalid") ;
        }
    }

    SLIP_FREE (option) ;
}

