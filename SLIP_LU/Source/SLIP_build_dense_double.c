//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_double: build dense double matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from double input */

#define SLIP_FREE_WORKSPACE                 \
    (*A_handle) = NULL ;                    \
    SLIP_delete_dense (&A) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_double
(
    SLIP_dense **A_handle,      // Dense matrix to construct
    double **b,                 // Set of values as doubles
    int32_t m,                  // number of rows
    int32_t n,                  // number of columns
    SLIP_options* option
)
{

    SLIP_info ok;
    if (!b || !A_handle)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }

    SLIP_dense *A = slip_create_dense ( ) ;
    if (A == NULL)
    {
        SLIP_FREE_WORKSPACE ;
        return (SLIP_OUT_OF_MEMORY) ;
    }

    SLIP_CHECK (slip_dense_alloc (A, m, n)) ;
    SLIP_CHECK (slip_expand_double_mat (A->x, b, A->scale, m, n, option)) ;

    (*A_handle) = A ;
    return (SLIP_OK) ;
}

