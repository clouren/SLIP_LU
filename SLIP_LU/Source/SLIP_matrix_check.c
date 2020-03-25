//------------------------------------------------------------------------------
// SLIP_LU/SLIP_matrix_check: check if a matrix is OK
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

/* check the validity of a SLIP_sparse sparse matrix in compressed-
   sparse column form.  Derived from SuiteSparse/MATLAB_TOOLS/spok. */

// option->print_level:
//      0: nothing
//      1: just errors
//      2: errors and terse output
//      3: verbose

SLIP_info SLIP_matrix_check     // returns a SLIP_LU status code
(
    SLIP_matrix *A,     // matrix to check
    SLIP_options* option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (option == NULL)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }
    int print_level = option->print_level;
    if (A == NULL)
    {
        if (print_level >= 1)
        {
            printf ("A is NULL.\n") ;
        }
        return (SLIP_INCORRECT_INPUT) ;
    }
    
    //TODO Should we change this to work on all matrix types?
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    int64_t i, j, p, pend ;

    int64_t m = A->m ;
    int64_t n = A->n ;
    int64_t nzmax = A->nzmax ;
    int64_t nz = A->nz ;            // also nz == Ap [n]
    int64_t *Ap = A->p ;
    int64_t *Ai = A->i ;
    mpz_t *Ax = A->x.mpz ;

    if (print_level >= 2)
    {
        printf ("SLIP_sparse: m %"PRId64" n %"PRId64" nz %"PRId64
            " nzmax %"PRId64"\n", m, n, nz, nzmax) ;
    }

    //--------------------------------------------------------------------------
    // check the dimensions
    //--------------------------------------------------------------------------

    if (m < 0)      // or (m <= 0) ?
    {
        if (print_level > 0) printf ("m invalid\n") ;
        return (SLIP_INCORRECT_INPUT) ;
    }
    if (n < 0)      // or (n <= 0) ?
    {
        if (print_level > 0) printf ("n invalid\n") ;
        return (SLIP_INCORRECT_INPUT) ;
    }
    if (nzmax < 0)      // or (nzmax <= 0) ?
    {
        if (print_level > 0) printf ("nzmax invalid\n") ;
        return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------
    // check the column pointers
    //--------------------------------------------------------------------------

    if (nz > 0 && (Ap == NULL || Ap [0] != 0 || Ap [n] != nz))
    {
        // column pointers invalid
        if (print_level > 0) printf ("p invalid\n") ;
        return (SLIP_INCORRECT_INPUT) ;
    }
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;
        pend = Ap [j+1] ;
        if (pend < p || pend > nz)
        {
            // column pointers not monotonically non-decreasing
            if (print_level > 0) printf ("p invalid\n") ;
            return (SLIP_INCORRECT_INPUT) ;
        }
    }

    //--------------------------------------------------------------------------
    // check the row indices and values
    //--------------------------------------------------------------------------

    if (nzmax > 0 && (Ai == NULL || Ax == NULL))
    {
        // row indices or values not present
        if (print_level > 0) printf ("i or x invalid\n") ;
        return (SLIP_INCORRECT_INPUT) ;
    }

    // allocate workspace to check for duplicates
    int64_t *mark = SLIP_calloc (n, sizeof (int64_t)) ;
    if (mark == NULL)
    {
        // out of memory
        if (print_level > 0) printf ("out of memory\n") ;
        return (SLIP_OUT_OF_MEMORY) ;
    }

    for (j = 0 ; j < n ; j++)
    {
        if (print_level >= 2)
        {
            printf ("column %"PRId64" :\n", j) ;
        }
        int64_t marked = j+1 ;
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            i = Ai [p] ;
            if (i < 0 || i >= m || mark [i] == marked)
            {
                // row indices out of range, or duplicate
                if (print_level > 0) printf ("invalid index\n") ;
                SLIP_FREE (mark) ;
                return (SLIP_INCORRECT_INPUT) ;
            }
            if (print_level >= 2)
            {
                printf ("  row %"PRId64" : ", i) ;
                SLIP_info status = SLIP_gmp_printf ( "%Zd " , Ax [p]) ;
                if (status < 0)
                {
                    SLIP_FREE (mark) ;
                    printf (" error: %d\n", status) ;
                    return (status) ;
                }
                printf ("\n") ;
            }
            mark [i] = marked ;
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE (mark) ;
    return (SLIP_OK) ;
}

