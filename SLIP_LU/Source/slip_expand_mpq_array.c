//------------------------------------------------------------------------------
// SLIP_LU/slip_expand_mpq_array: convert mpq array to mpz
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpq array of size n into an appropriate mpz
 * array of size n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq arrays to be used in SLIP LU 
 */

#define SLIP_FREE_WORKSPACE              \
    SLIP_delete_mpz_array(&x3, n);  \
    SLIP_delete_mpq_array(&x4, n);  \
    SLIP_MPZ_CLEAR(temp);

# include "SLIP_LU_internal.h"

SLIP_info slip_expand_mpq_array
(
    mpz_t* x_out,    //mpz array 
    mpq_t* x,     //mpq array that needs to be converted
    mpq_t scale,  //scaling factor. x_out = scale*x
    int32_t n     //size of x
)
{
    SLIP_info ok;
    mpz_t temp;
    SLIP_MPZ_SET_NULL(temp);
    ok = SLIP_mpz_init(temp);
    mpz_t* x3 = SLIP_create_mpz_array(n);    // Initialize arrays
    mpq_t* x4 = SLIP_create_mpq_array(n);
    if (!x3 || !x4 || ok != SLIP_OK) 
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    
    for (int32_t i = 0; i < n; i++)                  // x3 = denominators of x
    {
        SLIP_CHECK(SLIP_mpq_get_den(x3[i], x[i]));
    }
    SLIP_CHECK(SLIP_mpz_set(temp,x3[0]));
    
    // Find LCM of denominators of x
    for (int32_t i = 1; i < n; i++)
    {
        SLIP_CHECK(SLIP_mpz_lcm(temp, x3[i], temp));
    }
    SLIP_CHECK(SLIP_mpq_set_z(scale,temp));
    
    for (int32_t i = 0; i < n; i++)    
    {
        // x4[i] = x[i]*temp
        SLIP_CHECK(SLIP_mpq_mul(x4[i], x[i], scale));
    
        // x_out[i] = x4[i]
        SLIP_CHECK(SLIP_mpz_set_q(x_out[i], x4[i]));
    }
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

