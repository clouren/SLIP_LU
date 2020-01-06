//------------------------------------------------------------------------------
// SLIP_LU/slip_get_nonzero_pivot: find a nonzero pivot in a column
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* This function obtains the first eligible nonzero pivot
 * This is enabled if the user sets option->pivot = SLIP_FIRST_NONZERO
 * NOTE: This pivoting scheme is not recommended
 *
 * On output, the kth pivot is returned
 * 
 */
int32_t slip_get_nonzero_pivot // find the first eligible nonzero pivot
(
    int32_t *pivot,   // the index of first eligible nonzero pivot
    mpz_t* x,         // kth column of L and U
    int32_t* pivs,    // vector indicating which rows are pivotal
    int32_t n,        // size of x
    int32_t top,      // nonzero pattern is located in xi[top..n-1]
    int32_t* xi       // nonzero pattern of x
)
{
    if (!x || !pivs || !xi) {return SLIP_INCORRECT_INPUT;}
    int32_t inew, sgn;
    *pivot = -1; // used later to check for singular matrix
    SLIP_info ok;

    //--------------------------------------------------------------------------
    // Iterate accross the nonzeros in x
    //--------------------------------------------------------------------------
    for (int32_t i = top; i < n; i++)
    {
        // inew is the location of the ith nonzero
        inew = xi[i];
        // check if x[inew] is an eligible pivot
        ok = SLIP_mpz_sgn(&sgn, x[inew]);
        if (ok != SLIP_OK) {return ok;}
        if (sgn != 0 && pivs [inew] < 0)
        {
            *pivot = inew;
            // End the loop
            break;
        }
    }
    if (*pivot == -1)
    {
        return SLIP_SINGULAR;
    }
    else
    {
        return SLIP_OK;
    }
}
