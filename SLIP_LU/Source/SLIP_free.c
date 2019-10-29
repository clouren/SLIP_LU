//------------------------------------------------------------------------------
// SLIP_LU/SLIP_free: wrapper for free
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

/* Purpose: Free the memory associated with the pointer x
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

