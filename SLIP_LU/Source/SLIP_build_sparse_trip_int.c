//------------------------------------------------------------------------------
// SLIP_LU/SLIP_build_sparse_trip_int: build sparse matrix from int
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function will allow the user to take a matrix of their defined
 * type (in this case int) and convert it from their
 * triplet form to our data structure. The integrity of the user defined arrays
 * are maintained (therefore, one would need to delete these arrays)
 *
 * On output, the SLIP_sparse* A contains the user's matrix
 *
 */

#define SLIP_FREE_WORKSPACE                  \
    SLIP_delete_mpz_array(&x_new, nz);

 #include "SLIP_LU_internal.h"
 
 SLIP_info SLIP_build_sparse_trip_int
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    int32_t *x,         // Set of values in int
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x, I, and J vectors)
)
{
    SLIP_info ok;
    if (!I || !J || !A_output || !x || n <= 0 || nz <= 0 || !A_output->scale)
    {
        return SLIP_INCORRECT_INPUT;
    }

    mpz_t *x_new = SLIP_create_mpz_array(nz);
    if (x_new == NULL)
    {
        return SLIP_OUT_OF_MEMORY;
    }

    for (int32_t k = 0; k < nz; k++)
    {
        SLIP_CHECK(SLIP_mpz_set_si(x_new[k], x[k]));
    }
    SLIP_CHECK(SLIP_mpq_set_ui(A_output->scale, 1,1));

    SLIP_CHECK(slip_trip_to_mat(A_output, I, J, x_new, n, nz));

    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}
