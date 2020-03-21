//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_int64: build dense int64_t matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Purpose: Build a SLIP_dense matrix C from integer (int64_t) input 
// matrix B.
//
// The input is an m-by-n matrix B, held in (int64_t **) format.  That is, B
// is a pointer to an array of size m, containing (int64_t *) pointers.
// B[i] is a pointer to an array of size n, so that (i,j)th entry in B can
// be accessed as B[i][j].
//
// The output &C is a handle, or a pointer to the matrix C of size m-by-n.  The
// SLIP_dense matrix C is a pointer to a SLIP_dense struct, containing values
// C->x of type mpz_t, dimensions C->m and C->n, and a scale factor C->scale.

// TODO: why int64_t as a datatype?  We should use int64_t

#define SLIP_FREE_ALL                 \
    (*C_handle) = NULL ;                    \
    SLIP_delete_dense (&C) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_int64
(
    SLIP_dense **C_handle,      // Dense matrix to construct
    // inputs, not modified:
    int64_t **B,                // Set of values as ints
    int64_t m,                  // number of rows
    int64_t n                   // number of columns
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_copy (&C, SLIP_DENSE, SLIP_MPZ, B, option)
    // instead, and remove this function (or at least make non-user-callable).

    SLIP_info info ;
    if (!B || !C_handle)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------

    SLIP_dense *C = slip_create_dense ( ) ;
    if (C == NULL)
    {
        SLIP_FREE_ALL ;
        return (SLIP_OUT_OF_MEMORY) ;
    }

    SLIP_CHECK (slip_dense_alloc(C, m, n)) ;
    for (int64_t i = 0; i < m; i++)
    {
        for (int64_t j = 0; j < n; j++)
        {
            SLIP_CHECK(SLIP_mpz_set_si( C->x[i][j], B [i][j]));
        }
    }

    SLIP_CHECK(SLIP_mpq_set_ui(C->scale, 1, 1));

    (*C_handle) = C ;
    return (SLIP_OK) ;
}

