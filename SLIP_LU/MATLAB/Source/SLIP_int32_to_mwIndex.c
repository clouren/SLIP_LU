#include "SLIP_LU_mex.h"

/* Purpose: This function converts an int to an mwIndex */
void SLIP_int32_to_mwIndex
(
    mwIndex* y,    // the mwIndex vector
    int32_t* x,    // the int32_t vector
    int32_t n       // the size of x and y
)
{
    for (int32_t i = 0; i < n; i++)
    {        
        y[i] = (mwSize) x[i];
    }
}
