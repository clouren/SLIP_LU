//------------------------------------------------------------------------------
// SLIP_LU/slip_cast_matrix: create a dense typecasted matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// slip_cast_matrix constructs a dense nz-by-1 matrix Y that holds the
// typecasted values of the input matrix A.  The input matrix A can be of any
// kind (CSC, triplet, or dense) and any type.

#define SLIP_FREE_ALL                   \
    SLIP_matrix_free (&Y, option) ;

#include "SLIP_LU_internal.h"

SLIP_info slip_cast_matrix
(
    SLIP_matrix **Y_handle,     // nz-by-1 dense matrix to create
    SLIP_type Y_type,           // type of Y
    SLIP_matrix *A,             // matrix with nz entries
    SLIP_options *option        // Command options, if NULL defaults are used
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info = SLIP_OK ;
    // TODO NULL option always return nz=-1, but it's awkward to create default
    // option before basic input check
    int64_t nz = SLIP_matrix_nnz (A, option) ;
    SLIP_matrix *Y = NULL ;
    if (Y_handle == NULL || A == NULL ||  nz < 0)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }
    (*Y_handle) = NULL ;

    SLIP_options* option2 = NULL;
    if (option == NULL)
    {
        option2 = SLIP_create_default_options();
    }
    else
    {
        option2 = option;
    }
    
    //--------------------------------------------------------------------------
    // allocate Y (shallow if Y_type is the same as A->type)
    //--------------------------------------------------------------------------

    SLIP_CHECK (SLIP_matrix_allocate (&Y, SLIP_DENSE, Y_type,
        nz, 1, nz, Y_type == A->type, true, option2)) ;

    //--------------------------------------------------------------------------
    // typecast the values from A into Y
    //--------------------------------------------------------------------------

    if (Y_type == A->type)
    {

        //----------------------------------------------------------------------
        // Y is shallow; just copy in the pointer of the values of A 
        //----------------------------------------------------------------------

        switch (Y_type)
        {
            case SLIP_MPZ:   Y->x.mpz   = A->x.mpz   ;
            break;
            case SLIP_MPQ:   Y->x.mpq   = A->x.mpq   ;
            break;
            case SLIP_MPFR:  Y->x.mpfr  = A->x.mpfr  ;
            break;
            case SLIP_INT64: Y->x.int64 = A->x.int64 ;
            break;
            case SLIP_FP64:  Y->x.fp64  = A->x.fp64  ;
            break;
            default: SLIP_FREE_ALL ; return (SLIP_INCORRECT_INPUT) ;
        }

    }
    else
    {

        //----------------------------------------------------------------------
        // Y is deep; typecast the values from A into Y
        //----------------------------------------------------------------------

        SLIP_CHECK (slip_cast_array (SLIP_X (Y), Y->type,
            SLIP_X (A), A->type, nz, Y->scale, option2)) ;
           
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    /*TODO add this?
    if (option == NULL)
    {
        SLIP_FREE(option2);
    }
    */
    (*Y_handle) = Y;
    SLIP_CHECK (info) ;
    return (SLIP_OK) ;
}

