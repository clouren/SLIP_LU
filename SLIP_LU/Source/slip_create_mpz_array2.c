//------------------------------------------------------------------------------
// SLIP_LU/slip_create_mpz_array2: create an mpz array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function creates an mpz array of size n and allocates
 * memory for numbers of bit size prec. If the relative size of numbers is 
 * known ahead of time, this is more efficient than the
 * SLIP_create_mpz_array
 */
mpz_t* slip_create_mpz_array2
(
    int32_t n,     // size of the array
    int32_t size   // Relative size of numbers
)
{
    // Check input
    if (n <= 0 || size <= 0) {return NULL;}
    // Malloc space
    mpz_t* x = (mpz_t*) SLIP_calloc(n, SIZE_MPZ);
    if (!x) {return NULL;}
    for (int32_t i = 0; i < n; i++)
    {
        // Allocate x[i] for bit-length size
        if (SLIP_mpz_init2(x[i],size) != SLIP_OK)
        {
            // Out of memory
            SLIP_MPZ_SET_NULL(x[i]);
            SLIP_delete_mpz_array(&x, n);
            return NULL;
        }
    }
    return x;
}
