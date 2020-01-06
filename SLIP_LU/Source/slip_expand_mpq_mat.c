//------------------------------------------------------------------------------
// SLIP_LU/slip_expand_mpq_mat: convert mpq matrix to mpz
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts a mpq matrix of size m*n into an appropriate
 * mpz matrix of size m*n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq matrix to be used in SLIP LU 
 * 
 * on output, x_out is modified
 * 
 */

#define SLIP_FREE_WORKSPACE               \
    SLIP_delete_mpq_mat(&x4, m, n);  \
    SLIP_delete_mpz_mat(&x3, m, n);  \
    SLIP_MPZ_CLEAR(temp);

#include "SLIP_LU_internal.h"

SLIP_info slip_expand_mpq_mat
(
    mpz_t** x_out,// mpz mat
    mpq_t** x,    // mpq mat that needs to be converted
    mpq_t scale,  // scaling factor. x_out = scale*x
    int32_t m,    // number of rows of x
    int32_t n     // number of columns of x
)
{
    int32_t i, j;
    SLIP_info ok;
    mpq_t **x4 = NULL;
    mpz_t **x3 = NULL;
    mpz_t temp; SLIP_MPZ_SET_NULL(temp);
    x4 = SLIP_create_mpq_mat(m, n);
    x3 = SLIP_create_mpz_mat(m, n);
    ok = SLIP_mpz_init(temp);
    if (!x3 || !x4 || ok != SLIP_OK) 
    {
        // Out of memory
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // x3 = denominators of x
    //--------------------------------------------------------------------------
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            SLIP_CHECK(SLIP_mpq_get_den(x3[i][j], x[i][j]));
        }
    }
    SLIP_CHECK(SLIP_mpz_set(temp, x3[0][0]));
    
    //--------------------------------------------------------------------------
    // Find LCM of denominators of x
    //--------------------------------------------------------------------------
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            SLIP_CHECK(SLIP_mpz_lcm(temp, x3[i][j], temp));
        }
    }
    SLIP_CHECK(SLIP_mpq_set_z(scale, temp));
    
    //--------------------------------------------------------------------------
    // Make x integral
    //--------------------------------------------------------------------------
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            //x4[i][j] = x[i][j]*temp
            SLIP_CHECK(SLIP_mpq_mul(x4[i][j], x[i][j], scale));
        
            // x_out[i][j] = x4[i][j]
            SLIP_CHECK(SLIP_mpz_set_q(x_out[i][j], x4[i][j]));
        }
    }
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

