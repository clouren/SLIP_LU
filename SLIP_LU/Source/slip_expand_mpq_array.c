//------------------------------------------------------------------------------
// SLIP_LU/slip_expand_mpq_array: convert mpq array to mpz
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpq array of size n into an appropriate
 * mpz array of size n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq arrays to be used in SLIP LU.
 */

#define SLIP_FREE_ALL               \
    SLIP_MPZ_CLEAR(temp);           \
    SLIP_matrix_free(&x3, NULL);    \
    SLIP_matrix_free(&x4, NULL);    \

#include "SLIP_LU_internal.h"

SLIP_info slip_expand_mpq_array
(
    mpz_t* x_out,   // mpz array
    mpq_t* x,       // mpq array that needs to be converted
    mpq_t scale,    // scaling factor. x_out = scale*x
    int64_t n,      // size of x
    SLIP_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    
    if(!x_out || !x || n <0 || !option) return SLIP_INCORRECT_INPUT;
    
    SLIP_info info ;

    //--------------------------------------------------------------------------

    mpz_t temp;
    SLIP_matrix *x3 = NULL;
    SLIP_matrix *x4 = NULL;;
    SLIP_MPZ_SET_NULL(temp);
    SLIP_CHECK (SLIP_mpz_init(temp)) ;
        
    SLIP_matrix_allocate(&x3, SLIP_DENSE, SLIP_MPZ, n, 1, n, false, true, option);
    SLIP_matrix_allocate(&x4, SLIP_DENSE, SLIP_MPQ, n, 1, n, false, true, option);
    if (!x3 || !x4)
    {
        SLIP_FREE_ALL;
        return SLIP_OUT_OF_MEMORY;
    }

    for (int64_t i = 0; i < n; i++)                  // x3 = denominators of x
    {
        SLIP_CHECK(SLIP_mpq_get_den(x3->x.mpz[i], x[i]));
    }
    SLIP_CHECK(SLIP_mpz_set(temp,x3->x.mpz[0]));

    // Find LCM of denominators of x
    for (int64_t i = 1; i < n; i++)
    {
        SLIP_CHECK(SLIP_mpz_lcm(temp, x3->x.mpz[i], temp));
    }
    SLIP_CHECK(SLIP_mpq_set_z(scale,temp));

    for (int64_t i = 0; i < n; i++)
    {
        // x4[i] = x[i]*temp
        SLIP_CHECK(SLIP_mpq_mul(x4->x.mpq[i], x[i], scale));

        // x_out[i] = x4[i]
        SLIP_CHECK(SLIP_mpz_set_q(x_out[i], x4->x.mpq[i]));
    }
    SLIP_FREE_ALL;
    return SLIP_OK;
}

