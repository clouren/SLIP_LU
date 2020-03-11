//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_mpfr: build dense mpfr matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from mpfr_t input */

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_mpfr
(
    // TODO what does "allocated but unused yet" mean??
    SLIP_dense *A,          // Dense matrix, allocated but unused
    mpfr_t **b,             // Set of values as mpfr_t
    int32_t m,              // number of rows
    int32_t n,              // number of columns
    SLIP_options *option    // command options with precision for mpfr
)
{
    if (!b || !A ||!A->scale || !option)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A, m, n)) ;
    SLIP_CHECK(slip_expand_mpfr_mat(A->x, b, A->scale, m, n, option));
    return SLIP_OK;
}
