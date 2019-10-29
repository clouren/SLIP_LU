//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_int_mat: create dense int matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function creates an int matrix of size m*n. */
int32_t** SLIP_create_int_mat
(
    int32_t m,     // number of rows
    int32_t n      // number of columns
)
{
    // Check input
    if (m <= 0 || n <= 0) {return NULL;}
    int32_t **x = (int32_t**) SLIP_calloc(m, sizeof(int32_t*));
    if (!x) {return NULL;}
    for (int32_t i = 0; i < m; i++)     
    {
        x[i] = (int32_t*) SLIP_calloc(n, sizeof(int32_t));
        if (x[i] == NULL)
        {
            // Out of memory
            SLIP_delete_int_mat(&x, m, n);
            return NULL;
        }
    }
    return x;
}

