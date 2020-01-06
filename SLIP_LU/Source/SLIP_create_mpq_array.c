//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_mpq_array: create a dense mpq array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function creates an mpq array of size n. 
 * This function must be called for all mpq arrays created.
 */
mpq_t* SLIP_create_mpq_array
(
    int32_t n      // size of the array 
)        
{
    // Check input
    if (n <= 0) {return NULL;}
    // Malloc space
    mpq_t* x = (mpq_t*) SLIP_calloc(n, SIZE_MPQ);
    if (!x) {return NULL;}
    for (int32_t i = 0; i < n; i++)
    {
        if (SLIP_mpq_init(x[i]) != SLIP_OK)
        {
            // Out of memory
            SLIP_MPQ_SET_NULL(x[i]);
            SLIP_delete_mpq_array(&x,n);
            return NULL;
        }
    }
    return x;
}

