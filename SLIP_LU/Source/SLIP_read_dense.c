//------------------------------------------------------------------------------
// SLIP_LU/SLIP_read_dense: read a dense matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

/*
 * Purpose: This function reads in a RHS vector stored as a dense vector
 * This function requires the first line of the file as the number of row
 * and number of column, and the rest of the file lists each entry value.

 TODO: the file input format is a bad choice.  Use Matrix Market instead.

 */

SLIP_info SLIP_read_dense
(
    SLIP_dense *b,
    FILE* file          // file to read from (must already be open)
)
{
    if (b == NULL || file == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t nrows, ncols;
    SLIP_info ok;

    // Read in size of matrix & number of nonzeros
    ok = fscanf(file, "%d %d", &nrows, &ncols);
    if (feof(file) || ok < 2)
    {
        return SLIP_INCORRECT_INPUT;
    }
    if (slip_dense_alloc(b, nrows, ncols) != SLIP_OK)
    {
        return SLIP_OUT_OF_MEMORY;
    }

    for (int32_t i = 0; i < nrows; i++)
    {
        for (int32_t j = 0; j < ncols; j++)
        {
            ok = slip_gmp_fscanf(file, "%Zd", &(b->x[i][j]));
            if (ok < 1)
            {
                return SLIP_INCORRECT_INPUT;
            }
        }
    }
    return SLIP_OK;
}

