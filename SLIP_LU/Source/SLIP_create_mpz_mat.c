//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_mpz_mat: create a dense mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function creates a dense mpz_t matrix of size m*n to
 * default size
 */
mpz_t** SLIP_create_mpz_mat
(
    int32_t m,     // number of rows
    int32_t n      // number of columns
)
{
    // Check input
    if (m <= 0 || n <= 0) {return NULL;}
    // Calloc space
    mpz_t **x = (mpz_t**) SLIP_calloc(m, sizeof(mpz_t*));
    if (!x) {return NULL;}
    for (int32_t i = 0; i < m; i++)
    {
        x[i] = (mpz_t*) SLIP_calloc(n, SIZE_MPZ);
        if (x [i] == NULL)
        {
            // out of memory
            SLIP_delete_mpz_mat (&x, m, n) ;
            return (NULL) ;
        }
        for (int32_t j = 0; j < n; j++)
        {
            if (SLIP_mpz_init(x[i][j]) == SLIP_OUT_OF_MEMORY)
            {
                // out of memory
                SLIP_MPZ_SET_NULL(x[i][j]);
                SLIP_delete_mpz_mat (&x, m, n) ;
                return (NULL) ;
            }
        }
    }
    return x;
}

