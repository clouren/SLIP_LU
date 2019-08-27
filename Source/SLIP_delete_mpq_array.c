# include "SLIP_LU_internal.h"

/* Purpose: This function clears the memory used for an mpq vector of size n. 
 * Call this for all mpq vectors when done.
 * 
 * Input is a mpq_t** array which is destroyed upon function completion
 */
void SLIP_delete_mpq_array
(
    mpq_t** x,     // mpq array to be deleted
    int32_t n      // size of x 
)
{
    if (x == NULL || (*x) == NULL) {return;}
    for (int32_t i = 0; i < n; i++)
    {
        SLIP_MPQ_CLEAR((*x)[i]);
    }
    SLIP_FREE(*x);
} 
