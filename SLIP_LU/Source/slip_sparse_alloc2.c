//------------------------------------------------------------------------------
// SLIP_LU/slip_sparse_alloc2: allocate an uninitialized sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function allocates a SLIP LU matrix of size n*m with array
 * size nzmax. This version does not allocate individual the values in x. As a
 * result, it is more memory efficient, but also less user friendly.
 *
 * TODO: less user friendly?  It's not user-callable.  Rephrase this.
 *
 * See also slip_sparse_alloc.
 */

// TODO delete this function? use SLIP_matrix_allocate instead.
// TODO: does SLIP_matrix_allocate need a flag to tell it not to initialize
// the content of A->x.  This is only need by SLIP_LU_factorize, as an
// internal function, however.

// TODO: Perhaps the user-callable SLIP_matrix_allocate could simply be a
// 1-line function that calls slip_matrix_allocate2.  The user-accessible
// function pass in the flag to initialize all mpq, mpz, and mpfr entries.
// Then SLIP_LU_factorize could call slip_matrix_allocate2 directly, and pass
// in the flag to not initialize the content of these arrays.

#include "SLIP_LU_internal.h"

SLIP_info slip_sparse_alloc2
(
    SLIP_sparse* A, // sparse matrix data structure to be allocated
    int64_t n,      // number of columns
    int64_t m,      // number of rows (recall m=n assumed)
    int64_t nzmax   // size of allocated i and x arrays
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

    // Allocate memory for A->x, A->p, and A->i arrays
    A->x = (mpz_t*) SLIP_calloc(nzmax, sizeof(mpz_t));
    A->p = (int64_t*) SLIP_calloc(n+1, sizeof(int64_t));
    A->i = (int64_t*) SLIP_calloc(nzmax, sizeof(int64_t));

    if (!A->x || !A->p || !A->i) {return SLIP_OUT_OF_MEMORY;}
    return SLIP_OK;
}

