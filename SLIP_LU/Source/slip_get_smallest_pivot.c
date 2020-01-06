//------------------------------------------------------------------------------
// SLIP_LU/slip_get_smallest_pivot: find the smallest entry in a column
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function selects the pivot element as the smallest in the column 
 * This is activated by default or if the user sets
 * option->pivot = SLIP_SMALLEST
 * NOTE: This is the recommended pivoting scheme for SLIP LU
 * 
 * On output, the index of kth pivot is returned
 * 
 */

#define SLIP_FREE_WORKSPACE        \
    SLIP_MPZ_CLEAR(small);

# include "SLIP_LU_internal.h"

SLIP_info slip_get_smallest_pivot
(
    int32_t *pivot, // the index of smallest pivot
    mpz_t* x,       // kth column of L and U
    int32_t* pivs,  // vector indicating whether each row has been pivotal 
    int32_t n,      // dimension of problem
    int32_t top,    // nonzeros are stored in xi[top..n-1]
    int32_t* xi     // nonzero pattern of x
)
{
    // Check input
    if (!x || !pivs || !xi) {return SLIP_INCORRECT_INPUT;}
    // Flag is non-negative until we have an initial starting value for small
    int32_t i, inew, j, flag, sgn, r;
    *pivot = -1;
    SLIP_info ok;
    j = n;
    flag = top;
    mpz_t small; SLIP_MPZ_SET_NULL(small);
    SLIP_CHECK(SLIP_mpz_init(small));

    //--------------------------------------------------------------------------
    // Find an initial pivot. Fails if all terms are 0 in array x
    //--------------------------------------------------------------------------
    while (flag > -1 && flag < n)
    {
        // i location of first nonzero
        inew = xi[flag];
        
        //check if inew can be pivotal
        SLIP_CHECK(SLIP_mpz_sgn(&sgn, x[inew]));
        if (pivs[inew] < 0 && sgn != 0)
        {
            // Current smallest pivot
            SLIP_CHECK(SLIP_mpz_set(small, x[inew]));
            // Current smallest pivot location
            *pivot = inew;
            // Where to start the search for rest of nonzeros
            j = flag;
            // Exit the while loop
            flag = -5;
        }
        // Increment to next nonzero to search
        flag += 1;
    }

    //--------------------------------------------------------------------------
    // Iterate across remaining nonzeros
    //--------------------------------------------------------------------------
    for (i = j; i < n; i++)
    {
        inew = xi[i];
        // check if inew can be pivotal
        SLIP_CHECK(SLIP_mpz_cmpabs(&r, small, x[inew])); 
        if (pivs[inew] < 0 && r > 0)
        {
            SLIP_CHECK(SLIP_mpz_sgn(&sgn, x[inew]));
            if (sgn != 0)
            {
                // Current best pivot location
                *pivot = inew;
                // Current best pivot value
                SLIP_CHECK(SLIP_mpz_set(small, x[inew]));
            }
        }
    }
    // Clear memory for small
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

