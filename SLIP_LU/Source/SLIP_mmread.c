//------------------------------------------------------------------------------
// SLIP_LU/SLIP_mmread: read a sparse mpz matrix in Matrix Market format
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// TODO: this doesn't make any sense.  The file format here is NOT
// Matrix Market (where is the header??).

// TODO: move this to Demo folder and rename it?

#define SLIP_FREE_WORKSPACE           \
    SLIP_FREE(i);                     \
    SLIP_FREE(j);                     \
    SLIP_delete_mpz_array(&x_mpz, nz);

#include "SLIP_LU_internal.h"

/* Purpose: This function reads in a matrix stored in Matrix Market format */
SLIP_info SLIP_mmread
(
    SLIP_sparse* A,     // Matrix to be populated
    FILE* file          // file to read from (must already be open)
)
{
    SLIP_info ok;
    if (A == NULL || file == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    int32_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    ok = fscanf(file, "%d %d %d\n", &m, &n, &nz);
    if (feof(file) || ok < 3)
    {
        // TODO: untested
        return SLIP_INCORRECT_INPUT;
    }
    int32_t *i = (int32_t*) SLIP_malloc(nz * sizeof(int32_t));
    int32_t *j = (int32_t*) SLIP_malloc(nz * sizeof(int32_t));

    // Create an initialized input mpz vector
    mpz_t* x_mpz = SLIP_create_mpz_array(nz);
    if (!i || !j || !x_mpz)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    int32_t decrement;
    ok = slip_gmp_fscanf(file, "%d %d %Zd\n", &i[0], &j[0], &x_mpz[0]);
    if (feof(file) || ok < 3)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_INCORRECT_INPUT;
    }

    if (SLIP_MIN(i[0], j[0]) == 0)
    {
        // TODO this is invalid, and fragile; Matrix Market is always
        // one-based.  Where do you read in a file that is zero-based?
        decrement = 0;
    }
    else
    {
        decrement = 1;
        i[0]-=decrement;
        j[0]-=decrement;
    }

    // Read in the values from file
    for (int32_t p = 1; p < nz; p++)
    {
        ok = slip_gmp_fscanf(file, "%d %d %Zd\n", &i[p], &j[p], &x_mpz[p]);
	if ((feof(file) && p != nz-1) || ok < 3)
	{
	    SLIP_FREE_WORKSPACE;
	    return SLIP_INCORRECT_INPUT;
	}
        // Conversion from 1 based to 0 based
        i[p] -= decrement;
        j[p] -= decrement;
    }

    // Convert from triplet form to ccf
    ok = slip_trip_to_mat(A, i, j, x_mpz, n, nz);
    SLIP_FREE_WORKSPACE;
    return ok;
}

