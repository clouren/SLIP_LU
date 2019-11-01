//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_ccf_double: build sparse matrix from double
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

# include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------

#undef  SLIP_FREE_WORKSPACE
#define SLIP_FREE_WORKSPACE                 \
    SLIP_delete_mpz_array(&x_new, nz);
    
//------------------------------------------------------------------------------
// SLIP_build_sparse_ccf_mpz: build sparse matrix of mpz values
//------------------------------------------------------------------------------

/* Purpose: This function will allow the user to take a matrix of their defined
 * type (either double, mpfr_t, mpz_t, or mpq_t) and convert it from their
 * version of compressed column form to our data structure. The integrity of the
 * user defined arrays are maintained (therefore, one would need to delete these
 * arrays)
 *
 * On output, the SLIP_sparse* A structure contains the input matrix
 *
 */

SLIP_info SLIP_build_sparse_ccf_mpz
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    mpz_t *x,             // Set of values in full precision int.
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
)
{
    SLIP_info ok;
    if (!p || !I || !x || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    ok = slip_mpz_populate_mat(A_output, I, p, x, n, nz);
    if (ok != SLIP_OK) {return ok;}
    ok = slip_mpq_set_ui(A_output->scale, 1, 1);
    if (ok != SLIP_OK) {return ok;}
    return SLIP_OK;
}
