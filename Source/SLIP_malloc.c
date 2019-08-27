#include "SLIP_LU_internal.h"

/* Purpose: Define malloc and free for SLIP LU
 * 
 * Output arguments are not modified, returned is either a pointer to 
 * size space or a NULL pointer in the case of failure * 
 */
void* SLIP_malloc
(
    size_t size        // Size to alloc
)
{
    // ensure at least one byte is malloc'd
    if (size <= 0) {size = 1 ;}

    return (SLIP_MEMORY_MALLOC (size)) ;

}
