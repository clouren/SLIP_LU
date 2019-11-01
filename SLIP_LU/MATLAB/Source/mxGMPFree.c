#include "SLIP_LU_mex.h"

/* Purpose: A GMP free function. This allows GMP to use
 * MATLAB's mxFree instead of free */

// A GMP realloc function
void mxGMPFree 
(
    void* x,    // void* to be freed
    size_t a    // Size
)
{
    if (x) 
        mxFree(x);
}