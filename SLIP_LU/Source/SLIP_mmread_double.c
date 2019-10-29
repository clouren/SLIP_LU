#include "SLIP_LU_internal.h"


#define SLIP_FREE_WORKSPACE             \
    SLIP_FREE(i);                       \
    SLIP_FREE(j);                       \
    SLIP_FREE(x_doub);                  \
    SLIP_delete_mpz_array(&x_mpz, nz);


/*
 * Purpose: This function reads in a matrix stored in matrix market format
 * in this case the matrix is stored as a collection of doubles
 */
SLIP_info SLIP_mmread_double
(
    SLIP_sparse* A,        // Matrix to be populated
    FILE* file          // file to read from (must already be open)
)
{
    SLIP_info ok;
    if (A == NULL || file == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }
    // Read in triplet form first
    int32_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    ok = fscanf(file, "%d %d %d\n", &m, &n, &nz);
    if (feof(file) || ok < 3)
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    int32_t *i = (int32_t*) SLIP_malloc(nz* sizeof(int32_t));
    int32_t *j = (int32_t*) SLIP_malloc(nz* sizeof(int32_t));
    double *x_doub = (double*) SLIP_malloc(nz* sizeof(double));
    mpz_t *x_mpz = SLIP_create_mpz_array(nz);

    if (!i || !j || !x_doub || !x_mpz)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    int32_t decrement;
    ok = fscanf(file, "%d %d %lf\n", &(i[0]), &(j[0]), &(x_doub[0]));
    if (feof(file) || ok < 3)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_INCORRECT_INPUT;
    }

    if (SLIP_MIN(i[0], j[0]) == 0)
    {
        decrement = 0;
    }
    else
    {
        decrement = 1;
        i[0]-=decrement;
        j[0]-=decrement;
    }

    // Read in the values from file
    for (int32_t k = 1; k < nz; k++)
    {
        ok = fscanf(file, "%d %d %lf\n", &(i[k]), &(j[k]), &(x_doub[k]));
	if ((feof(file) && k != nz-1) || ok < 3)
	{
	    SLIP_FREE_WORKSPACE;
	    return SLIP_INCORRECT_INPUT;
	}
        // Conversion from 1 based to 0 based
        i[k] -= decrement;
        j[k] -= decrement;
    }

    // Convert x_doub from double to mpz_t
    ok = slip_expand_double_array(x_mpz, x_doub, A->scale, nz);
    // Convert from triplet form to ccf
    if (ok == SLIP_OK)
    {
        ok = slip_trip_to_mat(A, i, j, x_mpz, n, nz);
    }
    SLIP_FREE_WORKSPACE;
    return ok;
}

#undef SLIP_FREE_WORKSPACE
