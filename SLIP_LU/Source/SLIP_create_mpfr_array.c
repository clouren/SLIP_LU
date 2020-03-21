//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_mpfr_array: create a dense mpr array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates a MPFR array of desired precision. */

#include "SLIP_LU_internal.h"

mpfr_t* SLIP_create_mpfr_array
(
    int64_t n,     // size of the array
    SLIP_options *option// command options containing the prec for mpfr
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_allocate (&A, SLIP_DENSE, SLIP_MPFR, ...)
    // or make this an internal function.

    if (n <= 0) {return NULL;}

    //--------------------------------------------------------------------------

    mpfr_t* x = (mpfr_t*) SLIP_calloc(n, sizeof(mpfr_t));
    if (!x) {return NULL;}
    for (int64_t i = 0; i < n; i++)
    {
        if (SLIP_mpfr_init2(x[i], option->prec) != SLIP_OK)
        {
            SLIP_MPFR_SET_NULL(x[i]);
            SLIP_delete_mpfr_array(&x, n);
            return NULL;
        }
    }
    return x;
}

