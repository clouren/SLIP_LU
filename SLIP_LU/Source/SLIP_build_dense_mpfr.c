//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_mpfr: build dense mpfr matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// SLIP_build_dense_mpfr: Build a dense matrix from mpfr input
//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from mpfr_t input */
SLIP_info SLIP_build_dense_mpfr
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    mpfr_t **b,           // Set of values as mpfr_t
    int32_t m,            // number of rows
    int32_t n,            // number of columns
    SLIP_options *option  // command options containing the precision for mpfr
)
{
    if (!b || !A_output ||!A_output->scale || !option)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    SLIP_CHECK(slip_expand_mpfr_mat(A_output->x, b, A_output->scale,
        m, n, option));
    return SLIP_OK;
}
