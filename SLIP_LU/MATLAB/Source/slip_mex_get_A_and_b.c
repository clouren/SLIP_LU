//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_mex_get_A_and_b: Obtain user's A and b matrices
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"


/* Purpose: This function reads in the A matrix and right hand side vectors. */
void slip_mex_get_A_and_b
(
    SLIP_sparse *A,          // Internal SLIP Mat stored in ccf 
    SLIP_dense *b,           // mpz matrix used internally
    const mxArray* pargin[], // The input A matrix and options 
    int32_t nargin           // Number of input to the mexFunction
)
{
    if (!A || !pargin)
    {
        slip_mex_error (SLIP_INCORRECT_INPUT);
    }
    
    // Declare variables
    SLIP_info status;
    mxArray* tmp;
    mwSize nA, mA, nb, mb, Anz, k, j;
    mwIndex *Ap, *Ai;
    double *Ax, *bx;
    SLIP_options* option = SLIP_create_default_options();
    
    //--------------------------------------------------------------------------
    // Read in A
    //--------------------------------------------------------------------------
    // Read in Ap, Ai, Ax
    Ap =  mxGetJc(pargin[0]);
    Ai =  mxGetIr(pargin[0]);
    Ax =  mxGetDoubles(pargin[0]);
    if (!Ai || !Ap || !Ax)
    {
        slip_mex_error (SLIP_INCORRECT_INPUT);
    }

    // Get info about A
    nA = mxGetN (pargin[0]);
    mA = mxGetM (pargin[0]);
    Anz = (mwSize) Ap[nA];
    if (nA != mA)
    {
        mexErrMsgTxt("Error! A has to be square!");
    }

    int32_t* Ai_int = (int32_t*) SLIP_malloc(Anz* sizeof(int32_t));
    int32_t* Ap_int = (int32_t*) SLIP_malloc((nA+1)* sizeof(int32_t));
    if (!Ai_int || !Ap_int) {slip_mex_error(SLIP_OUT_OF_MEMORY);}

    // Does x have inf?
    slip_mex_check_for_inf(Ax, Anz);
    slip_mwIndex_to_int32(Ai_int, Ai, Anz);
    slip_mwIndex_to_int32(Ap_int, Ap, nA+1);

    // Input is already integral
    tmp = mxGetField(pargin[nargin-1], 0, "int");
    if (tmp == NULL)
    {
        mexErrMsgTxt("Error at getting the int parameter");
    }
    if (mxGetScalar(tmp) > 0) 
    {
        // Declare memory
        int32_t* Ax_int = (int32_t*) SLIP_malloc(Anz* sizeof(int32_t));
        if (!Ax_int) {slip_mex_error(SLIP_OUT_OF_MEMORY);}
        for (k = 0; k < Anz; k++)
        {
            Ax_int[k] = (int32_t) Ax[k];
        }
        // Create A with no scaling
        status = SLIP_build_sparse_ccf_int(A, Ap_int, Ai_int, Ax_int,
            (int32_t) nA, (int32_t) Anz);
        SLIP_FREE(Ax_int);
    }
    else
    {
        // Create A with scaling
        status = SLIP_build_sparse_ccf_double(A, Ap_int, Ai_int, Ax,
            (int32_t) nA, (int32_t) Anz, option);
    }
    if (status != SLIP_OK) 
    {
        mexErrMsgTxt("Issue reading in A. Please ensure matrix is correct "
            "and try again");
    }
    SLIP_FREE(Ai_int);
    SLIP_FREE(Ap_int);

    //--------------------------------------------------------------------------
    // Read in b
    //--------------------------------------------------------------------------
    if (nargin == 3)
    {
        bx =  mxGetDoubles(pargin[1]);
        if (!bx || !b)
        {
            slip_mex_error (SLIP_INCORRECT_INPUT);
        }
        // Get info about RHS vector(s)
        nb = mxGetN(pargin[1]);
        mb = mxGetM(pargin[1]);
        if (mb != mA)
        {
            mexErrMsgTxt("Error! Dimensions of A and B does not match!");
        }

        // Does b have inf?
        slip_mex_check_for_inf(bx, nb*mb);

        // Is b integral?
        tmp = mxGetField(pargin[nargin-1], 0, "intb");
        if (tmp == NULL)
        {   
            mexErrMsgTxt("Error at getting the intb parameter");
        }
        int32_t count = 0;
        if (mxGetScalar(tmp) > 0)
        {
            // populate bx to a int mat
            int32_t** bx_int = SLIP_create_int_mat((int32_t) mb, (int32_t) nb);
            if (!bx_int) {slip_mex_error(SLIP_OUT_OF_MEMORY);}
            for (k = 0; k < nb; k++)
            {
                for (j = 0; j < mb; j++)
                {
                    bx_int[j][k] = (int32_t) bx[count]; 
                    count++;
                }
            }

            // Create b
            status = SLIP_build_dense_int(b, bx_int, (int32_t) mb,
                (int32_t) nb);
            SLIP_delete_int_mat(&bx_int, (int32_t) mb, (int32_t) nb);
        }
        else
        {
            // populate bx to a double mat
            double** bx_doub;
            bx_doub = SLIP_create_double_mat((int32_t) mb, (int32_t) nb);
            if (!bx_doub) {slip_mex_error(SLIP_OUT_OF_MEMORY);}

            count = 0;
            for (k = 0; k < nb; k++)
            {
                for (j = 0; j < mb; j++)
                {
                    bx_doub[j][k] = bx[count];
                    count++;
                }
            }

            // Create b
            status = SLIP_build_dense_double(b, bx_doub, (int32_t) mb,
                (int32_t) nb, option);
            SLIP_delete_double_mat(&bx_doub, (int32_t) mb, (int32_t) nb);
        }
        if (status != SLIP_OK) 
        {
            mexErrMsgTxt("Issue reading in b. Please ensure RHS is correct "
                "and try again");
        }
    }
    SLIP_FREE(option);
}
