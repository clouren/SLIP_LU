//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_mpz: build dense mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from mpz input */

#define SLIP_FREE_WORKSPACE                 \
    (*A_handle) = NULL ;                    \
    SLIP_delete_dense (&A) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_mpz
(
    SLIP_dense **A_handle,      // Dense matrix to construct
    // inputs, not modified:
    mpz_t **b,                  // Set of values in full precision int.
    int32_t m,                  // number of rows
    int32_t n                   // number of columns
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

    SLIP_CHECK (slip_dense_alloc(A, m, n)) ;
    for (int32_t i = 0; i < m; i++)
    {
        for (int32_t j = 0; j < n; j++)
        {
            SLIP_CHECK (SLIP_mpz_set(A->x[i][j], b[i][j])) ;
        }
    }

    SLIP_CHECK (SLIP_mpq_set_ui (A->scale, 1, 1)) ;

    (*A_handle) = A ;
    return (SLIP_OK) ;
}

