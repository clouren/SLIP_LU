//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_mpfr_mat: create a dense mpfr matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function creates a mpfr_t matrix of size m*n with
 * precision prec
 */
mpfr_t** SLIP_create_mpfr_mat
(
    int32_t m,     // number of rows
    int32_t n,     // number of columns
    SLIP_options *option// command options containing the prec for mpfr
)
{
    // Check input
    if (m <= 0 || n <= 0) {return NULL;}

    mpfr_t **x = (mpfr_t**) SLIP_calloc(m, sizeof(mpfr_t*));
    if (!x) {return NULL;}
    for (int32_t i = 0; i < m; i++)
    {
        x[i] = NULL;
    }
    for (int32_t i = 0; i < m; i++)
    {
        x[i] = (mpfr_t*) SLIP_calloc(n, SIZE_MPFR);
        if (x[i] == NULL)
        {
            // Out of memory
            SLIP_delete_mpfr_mat(&x, m, n);
            return NULL;
        }
        for (int32_t j = 0; j < n; j++)
        {
            if (SLIP_mpfr_init2(x[i][j], option->prec) != SLIP_OK)
            {
                SLIP_MPFR_SET_NULL(x[i][j]);
                SLIP_delete_mpfr_mat(&x, m, n);
                return NULL;
            }
        }
    }
    return x;
}

