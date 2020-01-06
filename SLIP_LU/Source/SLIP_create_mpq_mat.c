//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_mpq_mat: create a dense mpq matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function creates a mpq_t matrix of size m*n. */
mpq_t** SLIP_create_mpq_mat
(
    int32_t m,     // number of rows
    int32_t n      // number of columns
)
{
    // Check input
    if (m <= 0 || n <= 0) {return NULL;}
    // Malloc space
    mpq_t **x = (mpq_t**) SLIP_calloc(m, sizeof(mpq_t*));
    if (!x) {return NULL;}
    for (int32_t i = 0; i < m; i++)
    {
        x[i] = (mpq_t*) SLIP_calloc(n, SIZE_MPQ);
        if ( x[i] == NULL)
        {
            // out of memory 
            SLIP_delete_mpq_mat(&x, m, n);
            return NULL;
        }
        for (int32_t j = 0; j < n; j++)
        {
            if (SLIP_mpq_init(x[i][j]) != SLIP_OK)
            {
                // Out of memory
                SLIP_MPQ_SET_NULL(x[i][j]);
                SLIP_delete_mpq_mat(&x, m, n);
                return NULL;
            }
        }            
    }
    return x;
}

