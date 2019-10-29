#include "SLIP_LU_internal.h"

/* Purpose: Free the memory associated with the pointer x
 * 
 * If we have defined MATLAB, then we use MATLAB's mxFree,
 * otherwise, we use default free.
 * 
 */


void SLIP_free
(
    void* x         // Pointer to be free'd
)
{
    if (x)
    {
        SLIP_MEMORY_FREE (x) ;
    }
}
