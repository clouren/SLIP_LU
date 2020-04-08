//------------------------------------------------------------------------------
// SLIP_LU/SLIP_calloc: wrapper for calloc
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Purpose: calloc space of size n*size.  Returns NULL on failure.

#include "slip_internal.h"

void* SLIP_calloc
(
    size_t n,          // Size of array
    size_t size        // Size of each entry
)
{
    // ensure at least one byte is calloc'd
    n = SLIP_MAX (n, 1) ;
    size = SLIP_MAX (size, 1) ;
    return (SLIP_MEMORY_CALLOC (n, size)) ;
}

