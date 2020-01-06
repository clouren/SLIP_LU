//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_mpz_array: create a dense mpz array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function creates an mpz array of size n and allocates
 * default size.
 */
mpz_t* SLIP_create_mpz_array 
(
    int32_t n      // Size of x
)        
{
    // Check inputs
    if (n <= 0) {return NULL;}
    
    // Malloc memory
    mpz_t* x = (mpz_t*) SLIP_calloc(n, SIZE_MPZ);
    if (!x) {return NULL;}
    for (int32_t i = 0; i < n; i++)
    {
        if (SLIP_mpz_init(x[i]) != SLIP_OK)
        {
            // Error, out of memory
            SLIP_MPZ_SET_NULL(x[i]);
            SLIP_delete_mpz_array(&x,n);
            return NULL;
        }
    }
    return x;
}

