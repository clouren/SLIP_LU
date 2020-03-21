//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_double_mat: create double dense matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function creates a double matrix of size m*n. */

#include "SLIP_LU_internal.h"

double** SLIP_create_double_mat
(
    int64_t m,     // number of rows (must be > 0)
    int64_t n      // number of columns (must be > 0)
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_allocate (&A, SLIP_DENSE, SLIP_DOUBLE, ...)

    if (m <= 0 || n <= 0) {return NULL;}

    //--------------------------------------------------------------------------

    double **x = (double**) SLIP_calloc(m, sizeof(double*));
    if (!x) {return NULL;}

    for (int64_t i = 0; i < m; i++)
    {
        x[i] = (double*) SLIP_calloc(n, sizeof(double));
        if (x[i] == NULL)
        {
            // Out of memory
            SLIP_delete_double_mat(&x, m, n);
            return NULL;
        }
    }
    return x;
}

