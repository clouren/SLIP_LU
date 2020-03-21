//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_mpz_array: create a dense mpz array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates an mpz array of size n and allocates default
 * size for each entry.
 */

#include "SLIP_LU_internal.h"

mpz_t* SLIP_create_mpz_array
(
    int64_t n      // Size of x
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_allocate (&A, SLIP_DENSE, SLIP_MPZ, ...)

    if (n <= 0) {return NULL;}

    //--------------------------------------------------------------------------

    // Malloc memory
    mpz_t* x = (mpz_t*) SLIP_calloc(n, sizeof(mpz_t));
    if (!x) {return NULL;}

    for (int64_t i = 0; i < n; i++)
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

