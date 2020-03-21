//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_mpfr: build dense mpfr matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Purpose: Build a SLIP_dense matrix C from mpfr_t input matrix B.
//
// The input is an m-by-n matrix B, held in (mpfr_t **) format.  That is, B
// is a pointer to an array of size m, containing (mpfr_t *) pointers.
// B[i] is a pointer to an array of size n, so that (i,j)th entry in B can
// be accessed as B[i][j].
//
// The output &C is a handle, or a pointer to the matrix C of size m-by-n.  The
// SLIP_dense matrix C is a pointer to a SLIP_dense struct, containing values
// C->x of type mpz_t, dimensions C->m and C->n, and a scale factor C->scale.

#define SLIP_FREE_ALL                  \
    (*C_handle) = NULL ;                    \
    SLIP_delete_dense (&C) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_mpfr
(
    SLIP_dense **C_handle,      // Dense matrix to construct
    // inputs, not modified:
    mpfr_t **B,             // Set of values as mpfr_t
    int64_t m,              // number of rows
    int64_t n,              // number of columns
    SLIP_options *option    // command options with precision for mpfr
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_copy (&C, SLIP_DENSE, SLIP_MPZ, B, option)
    // instead, and remove this function (or at least make non-user-callable).

    SLIP_info info ;
    if (!B || !C_handle || !option)
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

    SLIP_CHECK (slip_dense_alloc (C, m, n)) ;
    SLIP_CHECK (slip_expand_mpfr_mat (C->x, B, C->scale, m, n, option));

    (*C_handle) = C ;
    return (SLIP_OK) ;
}

