//------------------------------------------------------------------------------
// SLIP_LU/slip_sparse_alloc: allocate a sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

/* Purpose: This function allocates a SLIP LU matrix of size n*m with array
 * size nzmax. This function initializes each entry in A->x therefore they are
 * immediately ready to be operated on. This is less efficient but more user
 * friendly.
 *
 * TODO: more user friendly?  It's not user-callable.  Rephrase this.
 *
 * See also slip_sparse_alloc2.
 */

// TODO: delete this function; use SLIP_matrix_allocate instead

SLIP_info slip_sparse_alloc
(
    SLIP_sparse* A,// sparse matrix data structure to be allocated
    int32_t n,     // number of columns
    int32_t m,     // number of rows (recall m=n assumed)
    int32_t nzmax  // size of allocated i and x arrays
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (n <= 0 || m <= 0 || nzmax <= 0 || !A) {return SLIP_INCORRECT_INPUT;}

    //--------------------------------------------------------------------------

    A->m = m;                                    // Rows of the matrix
    A->n = n;                                    // Columns of the matrix
    A->nz = 0;                                   // Currently 0 nonzeros
    A->nzmax = nzmax;                            // Size of the vectors

    // Create and initialize A->x
    A->x = SLIP_create_mpz_array(nzmax);

    // Allocate A->p and A->i arrays.
    A->p = (int32_t*) SLIP_calloc(n+1, sizeof(int32_t));
    A->i = (int32_t*) SLIP_calloc(nzmax, sizeof(int32_t));

    if (!A->x || !A->p || !A->i) {return SLIP_OUT_OF_MEMORY;}
    return SLIP_OK;
}
