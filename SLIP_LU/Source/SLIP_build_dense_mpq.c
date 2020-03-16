//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_mpq: build dense mpq matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Purpose: Build a SLIP_dense matrix C from mpq_t input matrix B.
//
// The input is an m-by-n matrix B, held in (mpq_t **) format.  That is, B
// is a pointer to an array of size m, containing (mpq_t *) pointers.
// B[i] is a pointer to an array of size n, so that (i,j)th entry in B can
// be accessed as B[i][j].
//
// The output &C is a handle, or a pointer to the matrix C of size m-by-n.  The
// SLIP_dense matrix C is a pointer to a SLIP_dense struct, containing values
// C->x of type mpz_t, dimensions C->m and C->n, and a scale factor C->scale.

#define SLIP_FREE_WORKSPACE                 \
    (*C_handle) = NULL ;                    \
    SLIP_delete_dense (&C) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_mpq
(
    SLIP_dense **C_handle,      // Dense matrix to construct
    // inputs, not modified:
    mpq_t **B,                  // set of values as mpq_t
    int32_t m,                  // number of rows
    int32_t n                   // number of columns
)
{

    SLIP_info ok;
    if (!B || !C_handle)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }

    SLIP_dense *C = slip_create_dense ( ) ;
    if (C == NULL)
    {
        SLIP_FREE_WORKSPACE ;
        return (SLIP_OUT_OF_MEMORY) ;
    }

    SLIP_CHECK (slip_dense_alloc(C, m, n)) ;
    SLIP_CHECK (slip_expand_mpq_mat(C->x, B, C->scale, m, n));

    (*C_handle) = C ;
    return (SLIP_OK) ;
}

