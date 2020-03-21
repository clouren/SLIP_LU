//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_mpfr_array: delete an mpfr array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function clears the memory used for an mpfr array of size n.
 *
 * Input is a mpfr** array and its size. The input array is destroyed on
 * output.
 */

#include "SLIP_LU_internal.h"

void SLIP_delete_mpfr_array
(
    mpfr_t** x,    // mpfr array to be deleted
    int64_t n      // size of x
)
{

    // TODO: use SLIP_matrix_free (&A, option) ;
    // Move this functionality into SLIP_matrix_free.

    if (x == NULL || *x == NULL) {return;}
    for (int64_t i = 0; i < n; i++)
    {
        SLIP_MPFR_CLEAR((*x) [i]);
    }
    SLIP_FREE(*x);
}

