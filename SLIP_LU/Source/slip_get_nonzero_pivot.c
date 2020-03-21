//------------------------------------------------------------------------------
// SLIP_LU/slip_get_nonzero_pivot: find a nonzero pivot in a column
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function obtains the first eligible nonzero pivot
 * This is enabled if the user sets option->pivot = SLIP_FIRST_NONZERO
 *
 * NOTE: This pivoting scheme is NOT recommended for SLIP LU.  It is provided
 * for comparison with other pivoting options.
 *
 * On output, the kth pivot is returned.
 */

#include "SLIP_LU_internal.h"

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

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    if (!x || !pivs || !xi) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------

    int32_t inew, sgn;
    (*pivot) = -1; // used later to check for singular matrix

    //--------------------------------------------------------------------------
    // Iterate across the nonzeros in x
    //--------------------------------------------------------------------------

    for (int32_t i = top; i < n; i++)
    {
        // inew is the location of the ith nonzero
        inew = xi[i];
        // check if x[inew] is an eligible pivot
        SLIP_CHECK (SLIP_mpz_sgn (&sgn, x[inew])) ;
        if (sgn != 0 && pivs [inew] < 0)
        {
            (*pivot) = inew;
            // End the loop
            break;
        }
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    if ((*pivot) == -1)
    {
        return SLIP_SINGULAR;
    }
    else
    {
        return SLIP_OK;
    }
}

