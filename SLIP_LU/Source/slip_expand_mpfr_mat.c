//------------------------------------------------------------------------------
// SLIP_LU/slip_expand_mpfr_mat: convert mpfr matrix to mpz
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpfr matrix of size m*n and precision prec
 * to an appropriate mpz matrix of size m*n. To do this, the number is
 * multiplied by the appropriate power of 10 then the gcd is found. This
 * function allows mpfr arrays to be used within SLIP LU 
 */

#define SLIP_FREE_WORKSPACE          \
    SLIP_delete_mpfr_mat(&x3, m, n); \
    SLIP_MPFR_CLEAR(expon);          \
    SLIP_MPZ_CLEAR(temp_expon);      \
    SLIP_MPZ_CLEAR(gcd);             \
    SLIP_MPZ_CLEAR(one);             \
    SLIP_MPQ_CLEAR(temp);

# include "SLIP_LU_internal.h"

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
    int32_t i, j, k, l, r1, r2 = 1;
    bool nz_found = false;
    SLIP_info ok;
    mpfr_t expon, **x3 = NULL; SLIP_MPFR_SET_NULL(expon); 
    mpz_t temp_expon, gcd, one;
    SLIP_MPZ_SET_NULL(temp_expon);
    SLIP_MPZ_SET_NULL(gcd);
    SLIP_MPZ_SET_NULL(one);
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(SLIP_mpq_init(temp));
    SLIP_CHECK(SLIP_mpfr_init2(expon, option->prec));
    SLIP_CHECK(SLIP_mpz_init(temp_expon));
    SLIP_CHECK(SLIP_mpz_init(gcd));
    SLIP_CHECK(SLIP_mpz_init(one)); 
    x3 = SLIP_create_mpfr_mat(m, n, option);
    if (!x3)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
     
    // expon = 10^prec (overestimate)
    SLIP_CHECK(SLIP_mpfr_ui_pow_ui(expon, 10, option->prec, option->SLIP_MPFR_ROUND));
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            // x3[i][j] = x[i][j]*expon
            SLIP_CHECK(SLIP_mpfr_mul(x3[i][j], x[i][j], expon, option->SLIP_MPFR_ROUND));
            // x_out[i][j] = x3[i][j]
            SLIP_CHECK(SLIP_mpfr_get_z(x_out[i][j], x3[i][j], option->SLIP_MPFR_ROUND));
        }
    }
    
    SLIP_CHECK(SLIP_mpfr_get_z(temp_expon, expon, option->SLIP_MPFR_ROUND));
    SLIP_CHECK(SLIP_mpq_set_z(scale, temp_expon));
    
    //--------------------------------------------------------------------------
    // Find the gcd to reduce scale
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_mpz_set_ui(one, 1))
    // Find an initial GCD 
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (!nz_found)
            {
                SLIP_CHECK(SLIP_mpz_cmp_ui(&r1, x_out[i][j], 0));
                if (r1 != 0)
                {
                    nz_found = true;
                    k = i;
                    l = j;
                    SLIP_CHECK(SLIP_mpz_set(gcd, x_out[i][j]));;
                }
            }
            else
            {
                // Compute the GCD of the numbers
                SLIP_CHECK(SLIP_mpz_gcd(gcd, gcd, x_out[i][j]));
                SLIP_CHECK(SLIP_mpz_cmp(&r2, gcd, one));
                if (r2 == 0)
                {
                    break;
                }
            }
        }
        if (nz_found && r2 == 0)
        {
            break;
        }
    }
            
    if (!nz_found) // Entire matrix is zeros
    {
        SLIP_FREE_WORKSPACE;
        SLIP_mpq_set_z(scale, one);
        return SLIP_OK;
    }

    //--------------------------------------------------------------------------
    // Scale all entries to make as small as possible
    //--------------------------------------------------------------------------
    if (r2 != 0)         // If gcd == 1 stop
    {
        for (i = k; i < m; i++)
        {
            // start from the first found nonzero
            if (i == k) {j = l;}
            else        {j = 0;}
            for (; j < n; j++)
            {
                SLIP_CHECK(SLIP_mpz_divexact(x_out[i][j], x_out[i][j], gcd));
            }
        }
        SLIP_CHECK(SLIP_mpq_set_z(temp, gcd));
        SLIP_CHECK(SLIP_mpq_div(scale, scale, temp));
    }
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

