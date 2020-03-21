//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_int64_mat: create dense int64_t matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates an int64_t matrix of size m*n. */

#include "SLIP_LU_internal.h"

int64_t** SLIP_create_int64_mat
(
    int64_t m,     // number of rows (must be > 0)
    int64_t n      // number of columns (must be > 0)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_allocate (&A, SLIP_DENSE, SLIP_INT64, ...)

    if (m <= 0 || n <= 0) {return NULL;}

    //--------------------------------------------------------------------------

    int64_t **x = (int64_t**) SLIP_calloc(m, sizeof(int64_t*));
    if (!x) {return NULL;}

    for (int64_t i = 0; i < m; i++)
    {
        x[i] = (int64_t*) SLIP_calloc(n, sizeof(int64_t));
        if (x[i] == NULL)
        {
            // Out of memory
            SLIP_delete_int64_mat(&x, m, n);
            return NULL;
        }
    }
    return x;
}

