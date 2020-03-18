//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_int32: build dense int32_t matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Purpose: Build a SLIP_dense matrix C from integer (int32_t) input 
// matrix B.
//
// The input is an m-by-n matrix B, held in (int32_t **) format.  That is, B
// is a pointer to an array of size m, containing (int32_t *) pointers.
// B[i] is a pointer to an array of size n, so that (i,j)th entry in B can
// be accessed as B[i][j].
//
// The output &C is a handle, or a pointer to the matrix C of size m-by-n.  The
// SLIP_dense matrix C is a pointer to a SLIP_dense struct, containing values
// C->x of type mpz_t, dimensions C->m and C->n, and a scale factor C->scale.

#define SLIP_FREE_ALL                 \
    (*C_handle) = NULL ;                    \
    SLIP_delete_dense (&C) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_int32
(
    SLIP_dense **C_handle,      // Dense matrix to construct
    // inputs, not modified:
    int32_t **B,                // Set of values as ints
    int32_t m,                  // number of rows
    int32_t n                   // number of columns
)
{

    SLIP_info info ;
    if (!B || !C_handle)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }

    SLIP_dense *C = slip_create_dense ( ) ;
    if (C == NULL)
    {
        SLIP_FREE_ALL ;
        return (SLIP_OUT_OF_MEMORY) ;
    }

    SLIP_CHECK (slip_dense_alloc(C, m, n)) ;
    for (int32_t i = 0; i < m; i++)
    {
        for (int32_t j = 0; j < n; j++)
        {
            SLIP_CHECK(SLIP_mpz_set_si( C->x[i][j], B [i][j]));
        }
    }

    SLIP_CHECK(SLIP_mpq_set_ui(C->scale, 1, 1));

    (*C_handle) = C ;
    return (SLIP_OK) ;
}

