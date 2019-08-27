# include "SLIP_LU_internal.h"

/* 
 * Purpose: p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] in
 * to c 
 * This function is lightly modified from CSparse
 * 
 */
SLIP_info slip_cumsum
(
    int32_t *p,      // vector to store the sum of c
    int32_t *c,      // vector which is summed
    int32_t n        // size of c
)
{
    if (!p || !c || n <= 0) 
    {
        return SLIP_INCORRECT_INPUT;
    }

    int32_t i, nz = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        c [i] = p [i] ;
    }
    p [n] = nz ;
    return SLIP_OK ;
}
