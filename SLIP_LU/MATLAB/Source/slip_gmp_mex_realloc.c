#include "SLIP_LU_mex.h"

/* Purpose: A GMP reallocation function 
 * This allows GMP to use MATLAB's default realloc function 
 */

// A GMP realloc function
void* slip_gmp_mex_realloc 
(
    void* x,    // void* to be reallocated 
    size_t a,   // Previous size
    size_t b    // New size
)
{
    return mxRealloc(x,b);
}
