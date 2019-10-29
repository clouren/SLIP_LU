#include "SLIP_LU_mex.h"

/* Purpose: This function converts an mwIndex to an int32_t*/
void SLIP_mwIndex_to_int32
(
    int32_t* y,    // the int32_t vector
    mwIndex* x,    // the mwIndex vector
    mwSize n       // the size of x and y
)
{
    for (mwSize i = 0; i < n; i++)
    {  
        y[i] = (int32_t) x[i];
    }
}
