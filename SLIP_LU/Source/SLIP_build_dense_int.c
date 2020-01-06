//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_int: build dense int matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// SLIP_build_dense_int: Build a dense matrix from int input
//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from int input */
SLIP_info SLIP_build_dense_int
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    int32_t **b,          // Set of values as ints
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
    for (int32_t i = 0; i < m; i++)
    {
        for (int32_t j = 0; j < n; j++)
        {
            SLIP_CHECK(SLIP_mpz_set_si( A_output->x[i][j], b[i][j]));
        }
    }
    SLIP_CHECK(SLIP_mpq_set_ui(A_output->scale, 1, 1));
    return SLIP_OK;
}
