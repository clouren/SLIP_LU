//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_double: build dense double matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// SLIP_build_dense_double: Build a dense matrix from double input
//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from double input */
SLIP_info SLIP_build_dense_double
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    double **b,           // Set of values as doubles
    int32_t m,            // number of rows
    int32_t n,             // number of columns
    SLIP_options* option
)
{
    if (!b || !A_output || !A_output->scale)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    SLIP_CHECK(slip_expand_double_mat(A_output->x, b, A_output->scale, m, n, option));
    return SLIP_OK;
}
