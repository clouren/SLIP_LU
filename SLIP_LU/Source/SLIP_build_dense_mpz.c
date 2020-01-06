//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_mpz: build dense mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// SLIP_build_dense_mpz: Build a dense matrix from mpz input
//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from mpz input */
SLIP_info SLIP_build_dense_mpz
(
    SLIP_dense *A_output, // Dense matrix
    mpz_t **b,            // Set of values in full precision int.
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
            SLIP_CHECK(SLIP_mpz_set(A_output->x[i][j], b[i][j]));
        }
    }
    SLIP_CHECK(SLIP_mpq_set_ui(A_output->scale, 1, 1));
    return SLIP_OK;
}
