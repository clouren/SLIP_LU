//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_ccf_int: build sparse matrix from int
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#define SLIP_FREE_WORKSPACE                 \
    SLIP_delete_mpz_array(&x_new, nz);
    
# include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// SLIP_build_sparse_ccf_int: build sparse matrix from int values
//------------------------------------------------------------------------------

SLIP_info SLIP_build_sparse_ccf_int
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    int32_t *x,           // Set of values as doubles
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
)
{
    SLIP_info ok;
    if (!p || !I || !x || !A_output ||!A_output->scale)
    {
        return SLIP_INCORRECT_INPUT;
    }

    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (!x_new) {return SLIP_OUT_OF_MEMORY;}

    for (int32_t i = 0; i < nz; i++)
    {
            SLIP_CHECK(SLIP_mpz_set_si(x_new[i], x[i]));
    }
    SLIP_CHECK(SLIP_mpq_set_ui(A_output->scale, 1, 1));

    SLIP_CHECK(SLIP_mpz_populate_mat(A_output, I, p, x_new, n, nz));

    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}
