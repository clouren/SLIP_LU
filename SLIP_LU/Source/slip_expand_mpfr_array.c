//------------------------------------------------------------------------------
// SLIP_LU/slip_expand_mpfr_array: convert mprf aray to mpz
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpfr array of size n and precision prec to
 * an appropriate mpz array of size n. To do this, the number is multiplied by
 * the appropriate power of 10 then the gcd is found. This function allows mpfr
 * arrays to be used within SLIP LU 
 */

#define SLIP_FREE_WORKSPACE              \
    SLIP_MPFR_CLEAR(expon);         \
    SLIP_delete_mpfr_array(&x3, n); \
    SLIP_MPZ_CLEAR(temp_expon);     \
    SLIP_MPZ_CLEAR(gcd);            \
    SLIP_MPZ_CLEAR(one);            \
    SLIP_MPQ_CLEAR(temp);

# include "SLIP_LU_internal.h"

SLIP_info slip_expand_mpfr_array
(
    mpz_t* x_out,// full precision mpz array
    mpfr_t* x,  // mpfr array to be expanded
    mpq_t scale,// scaling factor used (x_out = scale*x)
    int32_t n,  // size of x
    SLIP_options *option  // command options containing the prec for mpfr
)
{   
    // Check input
    if (!x || !x_out || !scale || n <= 0) {return SLIP_INCORRECT_INPUT;}
    int32_t i, k, r1, r2 = 1;
    bool nz_found = false;
    SLIP_info ok;
    mpfr_t expon, *x3 = NULL; SLIP_MPFR_SET_NULL(expon);
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
    x3 = SLIP_create_mpfr_array(n, option);// Create the new x array
    if (!x3)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    // expon = 10^prec (overestimate)
    SLIP_CHECK(SLIP_mpfr_ui_pow_ui(expon, 10, option->prec, option->SLIP_MPFR_ROUND));
    for (i = 0; i < n; i++)
    {
        // x3[i] = x[i]*10^prec
        SLIP_CHECK(SLIP_mpfr_mul(x3[i], x[i], expon, option->SLIP_MPFR_ROUND));
        
        // x_out[i] = x3[i]
        SLIP_CHECK(SLIP_mpfr_get_z(x_out[i], x3[i], option->SLIP_MPFR_ROUND));
    }
    SLIP_CHECK(SLIP_mpfr_get_z(temp_expon, expon, option->SLIP_MPFR_ROUND));
    SLIP_CHECK(SLIP_mpq_set_z(scale, temp_expon));
    
    //--------------------------------------------------------------------------
    // Find the gcd to reduce scale 
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_mpz_set_ui(one, 1));
    // Find an initial GCD
    for (i = 0; i < n; i++)
    {
        if (!nz_found)
        {
            SLIP_CHECK(SLIP_mpz_cmp_ui(&r1, x_out[i], 0));
            if (r1 != 0)
            {
                nz_found = true;
                k = i;
                SLIP_CHECK(SLIP_mpz_set(gcd, x_out[i]));
            }
        }
        else
        {
            // Compute the GCD of the numbers, stop if gcd == 1
            SLIP_CHECK(SLIP_mpz_gcd(gcd, gcd, x_out[i]));
            SLIP_CHECK(SLIP_mpz_cmp(&r2, gcd, one));
            if (r2 == 0)
            {
                break;
            }
        }
    }

    if (!nz_found)     // Array is all zeros
    {
        SLIP_FREE_WORKSPACE;
        SLIP_mpq_set_z(scale, one);
        return SLIP_OK;
    }
    
    //--------------------------------------------------------------------------
    // Scale all entries to make as small as possible 
    //--------------------------------------------------------------------------
    if (r2 != 0)  // If gcd == 1 stop
    {
        for (i = k; i < n; i++)
        {
            SLIP_CHECK(SLIP_mpz_divexact(x_out[i],x_out[i],gcd));
        }
        SLIP_CHECK(SLIP_mpq_set_z(temp,gcd));
        SLIP_CHECK(SLIP_mpq_div(scale,scale,temp));
    }
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

