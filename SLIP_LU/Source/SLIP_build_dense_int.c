//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense_int: build dense int matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from int input */

#include "SLIP_LU_internal.h"

SLIP_info SLIP_build_dense_int
(
    // TODO what does "allocated but unused" mean??
    SLIP_dense *A,      // Dense matrix, allocated but unused
    int32_t **b,        // Set of values as ints
    int32_t m,          // number of rows
    int32_t n           // number of columns
)
{
    if (!b || !A ||!A->scale)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A, m, n)) ;
    for (int32_t i = 0; i < m; i++)
    {
        for (int32_t j = 0; j < n; j++)
        {
            SLIP_CHECK(SLIP_mpz_set_si( A->x[i][j], b[i][j]));
        }
    }
    SLIP_CHECK(SLIP_mpq_set_ui(A->scale, 1, 1));
    return SLIP_OK;
}
