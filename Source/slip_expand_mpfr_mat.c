# include "SLIP_LU_internal.h"


/* Purpose: This function converts a mpfr matrix of size m*n and precision prec
 * to an appropriate mpz matrix of size m*n. To do this, the number is
 * multiplied by the appropriate power of 10 then the gcd is found. This
 * function allows mpfr arrays to be used within SLIP LU 
 * NOTE: First element of input mpfr_t matrix must be nonzero
 */


#define SLIP_FREE_WORKSPACE               \
    SLIP_delete_mpfr_mat(&x3, m, n); \
    SLIP_MPFR_CLEAR(expon);          \
    SLIP_MPZ_CLEAR(temp_expon);      \
    SLIP_MPZ_CLEAR(gcd);             \
    SLIP_MPZ_CLEAR(one);             \
    SLIP_MPQ_CLEAR(temp);

SLIP_info slip_expand_mpfr_mat
(
    mpz_t** x_out,// mpz mat
    mpfr_t** x,   // mpfr matrix to be expanded
    mpq_t scale,  // scaling factor. x_out = scale*x
    int32_t m,    // number of rows of x
    int32_t n,    // number of columns of x
    SLIP_options *option// command options containing the prec for mpfr
)
{
    if (!x || !x_out || m <= 0 || n <= 0)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t i, j, r;
    SLIP_info ok;
    mpfr_t expon, **x3 = NULL; SLIP_MPFR_SET_NULL(expon); 
    mpz_t temp_expon, gcd, one;
    SLIP_MPZ_SET_NULL(temp_expon);
    SLIP_MPZ_SET_NULL(gcd);
    SLIP_MPZ_SET_NULL(one);
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(slip_mpq_init(temp));
    SLIP_CHECK(slip_mpfr_init2(expon, option->prec));
    SLIP_CHECK(slip_mpz_init(temp_expon));
    SLIP_CHECK(slip_mpz_init(gcd));
    SLIP_CHECK(slip_mpz_init(one)); 
    x3 = SLIP_create_mpfr_mat(m, n, option);
    if (!x3)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
     
    // expon = 10^prec (overestimate)
    SLIP_CHECK(slip_mpfr_ui_pow_ui(expon, 10, option->prec, MPFR_RNDN));
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            // x3[i][j] = x[i][j]*expon
            SLIP_CHECK(slip_mpfr_mul(x3[i][j], x[i][j], expon, MPFR_RNDN));
            // x_out[i][j] = x3[i][j]
            SLIP_CHECK(slip_mpfr_get_z(x_out[i][j], x3[i][j], MPFR_RNDN));
        }
    }
    
    SLIP_CHECK(slip_mpfr_get_z(temp_expon, expon, MPFR_RNDN));
    SLIP_CHECK(slip_mpq_set_z(scale, temp_expon));
    
    //--------------------------------------------------------------------------
    // Find the gcd to reduce scale
    //--------------------------------------------------------------------------
    // quit if x_out[0][0] == 0 (x is considered as incorrect input)
    SLIP_CHECK(slip_mpz_cmp_ui(&r, x_out[0][0], 0));
    if (r == 0)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_INCORRECT_INPUT;
    }
    SLIP_CHECK(slip_mpz_set(gcd, x_out[0][0]));
    SLIP_CHECK(slip_mpz_set_ui(one, 1));
    // Compute the GCD of the numbers
    for (i = 0; i < m && r != 0; i++)
    {
        for (j = 0; j < n && r != 0; j++)
        {
            SLIP_CHECK(slip_mpz_gcd(gcd, gcd, x_out[i][j]));
            SLIP_CHECK(slip_mpz_cmp(&r, gcd, one));
        }
    }
    
    //--------------------------------------------------------------------------
    // Scale all entries to make as small as possible
    //--------------------------------------------------------------------------
    if (r != 0)         // If gcd == 1 stop
    {
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                SLIP_CHECK(slip_mpz_divexact(x_out[i][j], x_out[i][j], gcd));
            }
        }
        SLIP_CHECK(slip_mpq_set_z(temp, gcd));
        SLIP_CHECK(slip_mpq_div(scale, scale, temp));
    }
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}
#undef SLIP_FREE_WORKSPACE
