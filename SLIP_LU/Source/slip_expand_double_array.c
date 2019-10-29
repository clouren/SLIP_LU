# include "SLIP_LU_internal.h"

/* Purpose: This function converts a double array of size n to an appropriate
 * mpz array of size n. To do this, the number is multiplied by 10^17 then, the
 * GCD is found. This function allows the use of matrices in double precision to
 * work with SLIP LU
 * NOTE: First element of input double array must be nonzero
 */

#define SLIP_FREE_WORKSPACE              \
    SLIP_delete_mpfr_array(&x3, n); \
    SLIP_MPZ_CLEAR(gcd);            \
    SLIP_MPZ_CLEAR(one);            \
    SLIP_MPQ_CLEAR(temp);



SLIP_info slip_expand_double_array
(
    mpz_t* x_out,//integral final array
    double* x,  //double array that needs to be made integral
    mpq_t scale,//the scaling factor used (x_out = scale * x)
    int32_t n   //size of x
)
{
    if (!x || !x_out || !scale || n <= 0)
    {
        return SLIP_INCORRECT_INPUT;
    }
    // int to store the comparison result
    int32_t r;
    SLIP_info ok;
    // Double precision accurate ~17 decimals
    double expon = pow(10, 17);
    // Quad precision in case input is huge
    mpfr_t* x3 = NULL;
    mpz_t gcd, one; SLIP_MPZ_SET_NULL(gcd); SLIP_MPZ_SET_NULL(one);
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(slip_mpq_init(temp));
    SLIP_CHECK(slip_mpz_init(gcd));
    SLIP_CHECK(slip_mpz_init(one));
    SLIP_options *option = SLIP_create_default_options();
    if (!option)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    x3 = SLIP_create_mpfr_array(n, option);
    SLIP_FREE(option);
    if (!x3)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    SLIP_CHECK(slip_mpq_set_d(scale, expon));           // scale = 10^17
    for (int32_t i = 0; i < n; i++)
    {
        // Set x3[i] = x[i]
        SLIP_CHECK(slip_mpfr_set_d(x3[i], x[i], MPFR_RNDN));

        // x3[i] = x[i] * 10^17
        SLIP_CHECK(slip_mpfr_mul_d(x3[i], x3[i], expon, MPFR_RNDN));

        // x_out[i] = x3[i]
        SLIP_CHECK(slip_mpfr_get_z(x_out[i], x3[i], MPFR_RNDN));
    }

    //--------------------------------------------------------------------------
    // Compute the GCD to reduce the size of scale
    //--------------------------------------------------------------------------
    // quit if x_out[0] == 0 (x is considered as incorrect input)
    SLIP_CHECK(slip_mpz_cmp_ui(&r, x_out[0], 0));
    if (r == 0)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_INCORRECT_INPUT;
    }
    SLIP_CHECK(slip_mpz_set(gcd, x_out[0]));
    SLIP_CHECK(slip_mpz_set_ui(one, 1));
    // Compute the GCD, stop if gcd == 1
    for (int32_t i = 1; i < n && r != 0; i++)
    {
        SLIP_CHECK(slip_mpz_gcd(gcd, gcd, x_out[i]));
        SLIP_CHECK(slip_mpz_cmp(&r, gcd, one));
    }

    //--------------------------------------------------------------------------
    // Scale all entries to make as small as possible
    //--------------------------------------------------------------------------
    if (r != 0 )             // If gcd == 1 then stop
    {
        for (int32_t i = 0; i < n; i++)
        {
            SLIP_CHECK(slip_mpz_divexact(x_out[i], x_out[i], gcd));
        }
        SLIP_CHECK(slip_mpq_set_z(temp, gcd));
        SLIP_CHECK(slip_mpq_div(scale, scale, temp));
    }
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}
#undef SLIP_FREE_WORKSPACE
