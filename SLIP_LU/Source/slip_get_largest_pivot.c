//------------------------------------------------------------------------------
// SLIP_LU/slip_get_largest_pivot: find a pivot entry in a column
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function selects the pivot element as the largest in the column
 * This is activated if the user sets option->pivot = SLIP_LARGEST
 * NOTE: This pivoting scheme is NOT recommended for SLIP LU
 * 
 * On output, the index of the largest pivot is returned
 * 
 */

#define SLIP_FREE_WORKSPACE   \
    SLIP_MPZ_CLEAR(big);

#include "SLIP_LU_internal.h"

SLIP_info slip_get_largest_pivot 
(
    int32_t *pivot,  // the index of largest pivot
    mpz_t* x,        // kth column of L and U
    int32_t* pivs,   // vector which indicates whether each row has been pivotal
    int32_t n,       // dimension of problem
    int32_t top,     // nonzero pattern is located in xi[top..n-1] 
    int32_t* xi      // nonzero pattern of x
)
{
    if (!x || !pivs || !xi) {return SLIP_INCORRECT_INPUT;}
    int32_t i, inew, r;
    SLIP_info ok;
    *pivot = -1;
    mpz_t big; SLIP_MPZ_SET_NULL(big);
    SLIP_CHECK(SLIP_mpz_init(big));

    //--------------------------------------------------------------------------
    // Iterate accross the nonzeros in x
    //--------------------------------------------------------------------------
    for (i = top; i < n; i++)
    {
        // Location of the ith nonzero
        inew = xi[i];
        // inew can be pivotal
        SLIP_CHECK(SLIP_mpz_cmpabs(&r, big, x[inew]));
        if (pivs[inew] < 0 && r < 0)
        {
            // Current largest pivot location
            *pivot = inew;
            // Current largest pivot value
            SLIP_CHECK(SLIP_mpz_set(big, x[inew]));
        }
    }

    // Frees the memory occupied by the pivot value
    SLIP_FREE_WORKSPACE;
    if (*pivot == -1)
    {
        return SLIP_SINGULAR;
    }
    else
    {
        return SLIP_OK;
    }
}

