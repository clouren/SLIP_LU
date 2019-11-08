//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_mpq: build dense mpq matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// SLIP_build_dense_mpq: Build a dense matrix from mpq input
//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from mpq_t input */
SLIP_info SLIP_build_dense_mpq
(
    SLIP_dense *A_output, // dense matrix, allocated but unused
    mpq_t **b,            // set of values as mpq_t
    int32_t m,            // number of rows
    int32_t n             // number of columns
)
{
    if (!b || !A_output ||!A_output->scale)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    SLIP_CHECK(slip_expand_mpq_mat(A_output->x, b, A_output->scale, m, n));
    return SLIP_OK;
}
