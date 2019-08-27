# include "SLIP_LU_internal.h"

/* 
 * Purpose: calloc space of size n*size
 * on failure, NULL is returned
 */

void* SLIP_calloc
(
    size_t n,          // Size of array
    size_t size        // Size to alloc
)
{
    // ensure at least one byte is calloc'd
    if (n <= 0) {n = 1 ;}
    if (size <= 0) {size = 1 ;}
    return (SLIP_MEMORY_CALLOC (n, size)) ;
}
