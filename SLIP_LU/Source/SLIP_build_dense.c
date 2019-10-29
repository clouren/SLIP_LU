//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_dense: build dense mpz, mpq, int, mprf, or double matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

// TODO split into 5 files?

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
    if (!b || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    for (int32_t i = 0; i < m; i++)
    {
        for (int32_t j = 0; j < n; j++)
        {
            SLIP_CHECK(slip_mpz_set(A_output->x[i][j], b[i][j]));
        }
    }
    SLIP_CHECK(slip_mpq_set_ui(A_output->scale, 1, 1));
    return SLIP_OK;
}

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
    if (!b || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    SLIP_CHECK(slip_expand_mpq_mat(A_output->x, b, A_output->scale, m, n));
    return SLIP_OK;
}

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
    if (!b || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    for (int32_t i = 0; i < m; i++)
    {
        for (int32_t j = 0; j < n; j++)
        {
            SLIP_CHECK(slip_mpz_set_si( A_output->x[i][j], b[i][j]));
        }
    }
    SLIP_CHECK(slip_mpq_set_ui(A_output->scale, 1, 1));
    return SLIP_OK;
}

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
    SLIP_options *option  // command options containing the prec for mpfr
)
{
    if (!b || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    SLIP_CHECK(slip_expand_mpfr_mat(A_output->x, b, A_output->scale,
        m, n, option));
    return SLIP_OK;
}

//------------------------------------------------------------------------------
// SLIP_build_dense_double: Build a dense matrix from double input
//------------------------------------------------------------------------------

/* Purpose: Build a dense matrix from double input */
SLIP_info SLIP_build_dense_double
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    double **b,           // Set of values as doubles
    int32_t m,            // number of rows
    int32_t n             // number of columns
)
{
    if (!b || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    SLIP_info ok;
    SLIP_CHECK (slip_dense_alloc(A_output, m, n)) ;
    SLIP_CHECK(slip_expand_double_mat(A_output->x, b, A_output->scale, m, n));
    return SLIP_OK;
}

