//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_mpz_array: delete a dense mpz array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function clears the memory used for an mpz vector of size n. 
 * Call this function for all mpz vectors when done.
 * 
 * Input is a mpz_t** array which is destroyed upon function completion
 * 
 */
void SLIP_delete_mpz_array
(
    mpz_t **x,      // mpz array to be deleted
    int32_t n       // Size of x
)
{
    if (x == NULL || (*x) == NULL) {return;}
    for (int32_t i = 0; i < n; i++)
    {
        if ( (*x) [i] != NULL)
        {
            SLIP_MPZ_CLEAR((*x) [i]);
        }
    }
    SLIP_FREE ((*x));
}

