//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_ccf_double: build sparse matrix from double
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#define SLIP_FREE_WORKSPACE                 \
    SLIP_delete_mpz_array(&x_new, nz);

# include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// SLIP_build_sparse_ccf_double: build sparse matrix from double values
//------------------------------------------------------------------------------

SLIP_info SLIP_build_sparse_ccf_double
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    double *x,            // Set of values as doubles
    int32_t n,            // dimension of the matrix
    int32_t nz,            // number of nonzeros in A (size of x and I vectors)
    SLIP_options* option
)
{
    SLIP_info ok;
    if (!p || !I || !x || !A_output || !A_output->scale)
    {
        return SLIP_INCORRECT_INPUT;
    }

    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (!x_new) {return SLIP_OUT_OF_MEMORY;}

    SLIP_CHECK(slip_expand_double_array(x_new, x, A_output->scale, nz, option));

    // Create our matrix
    SLIP_CHECK(SLIP_mpz_populate_mat(A_output, I, p, x_new, n, nz));

    // Free memory
    SLIP_FREE_WORKSPACE;
    // Success
    return SLIP_OK;
}
