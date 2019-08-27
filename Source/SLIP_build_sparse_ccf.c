# include "SLIP_LU_internal.h"

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

#define SLIP_FREE_WORKSPACE                 \
    SLIP_delete_mpz_array(&x_new, nz);

SLIP_info SLIP_build_sparse_ccf_double
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    double *x,            // Set of values as doubles
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
)
{
    SLIP_info ok;
    if (!p || !I || !x || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (!x_new) {return SLIP_OUT_OF_MEMORY;}

    SLIP_CHECK(slip_expand_double_array(x_new, x, A_output->scale, nz));

    // Create our matrix
    SLIP_CHECK(slip_mpz_populate_mat(A_output, I, p, x_new, n, nz));

    // Free memory
    SLIP_FREE_WORKSPACE;
    // Success
    return SLIP_OK;
}

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
    if (!p || !I || !x || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (!x_new) {return SLIP_OUT_OF_MEMORY;}

    for (int32_t i = 0; i < nz; i++)
    {
            SLIP_CHECK(slip_mpz_set_si(x_new[i], x[i]));
    }
    SLIP_CHECK(slip_mpq_set_ui(A_output->scale, 1, 1));

    SLIP_CHECK(slip_mpz_populate_mat(A_output, I, p, x_new, n, nz));

    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

SLIP_info SLIP_build_sparse_ccf_mpq
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    mpq_t *x,             // Set of values as mpq_t rational numbers
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
)
{
    SLIP_info ok;
    if (!p || !I || !x || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (!x_new) {return SLIP_OUT_OF_MEMORY;}

    SLIP_CHECK(slip_expand_mpq_array(x_new, x, A_output->scale, nz));

    SLIP_CHECK(slip_mpz_populate_mat(A_output, I, p, x_new, n, nz));

    SLIP_FREE_WORKSPACE;

    return SLIP_OK;
}

SLIP_info SLIP_build_sparse_ccf_mpfr
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    mpfr_t *x,            // Set of values as doubles
    int32_t n,            // dimension of the matrix
    int32_t nz,           // number of nonzeros in A (size of x and I vectors)
    SLIP_options *option  // command options containing the prec for mpfr
)
{
    SLIP_info ok;
    if (!p || !I || !x || !A_output)
    {
        return SLIP_INCORRECT_INPUT;
    }

    mpz_t* x_new = SLIP_create_mpz_array(nz);
    if (!x_new) {return SLIP_OUT_OF_MEMORY;}

    SLIP_CHECK(slip_expand_mpfr_array(x_new, x, A_output->scale, nz, option));

    SLIP_CHECK(slip_mpz_populate_mat(A_output, I, p, x_new, n, nz));

    SLIP_FREE_WORKSPACE;

    return SLIP_OK;
}

#undef SLIP_FREE_WORKSPACE
