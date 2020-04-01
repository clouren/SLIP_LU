//------------------------------------------------------------------------------
// SLIP_LU/SLIP_matrix_check: check if a matrix is OK
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------
#define SLIP_FREE_ALL    \
    SLIP_FREE(mark);     \
    SLIP_FREE(bmark);

#include "slip_LU_internal.h"

/* check the validity of a SLIP_matrix */

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

    int print_level = SLIP_GET_PRINT_LEVEL(option);
    if (A == NULL)
    {
        if (print_level > 0)
        {
            printf ("A is NULL.\n") ;
        }
        return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------
    // get the input matrix
    //--------------------------------------------------------------------------

    int64_t i, j, p, pend ;
    int64_t* mark = NULL;   // used for checking duplicate index of CSC
    bool* bmark = NULL;     // used for checking duplicate index of triplet

    int64_t m = A->m ;
    int64_t n = A->n ;
    int64_t nzmax = A->nzmax ;
    int64_t nz = SLIP_matrix_nnz(A, option);    // Number of nonzeros in A

    if (A->kind < 0 || A->kind > 2 || A->type < 0 || A->type > 4)
    {
        if (print_level > 0)
        {
            printf ("A has invalid kind or type.\n") ;
        }
        return (SLIP_INCORRECT_INPUT) ;
    }
    else
    {
        if (print_level > 1)
        {
            printf ("SLIP_matrix: nrows: %"PRId64", ncols: %"PRId64", nz:"
            "%"PRId64", nzmax: %"PRId64", kind: %s, type: %s\n", m, n, nz,
            nzmax, A->kind < 1 ? "CSC" : A->kind < 2 ? "Triplet" : "Dense",
            A->type < 1 ? "MPZ" : A->type < 2 ? "MPQ" : A->type < 3 ?
            "MPFR" : A->type < 4 ? "int64" : "double") ;
        }
    }

    //--------------------------------------------------------------------------
    // check the dimensions
    //--------------------------------------------------------------------------

    if (m < 0)      // or (m <= 0) ?
    {
        if (print_level > 0) {printf ("m invalid\n") ;}
        return (SLIP_INCORRECT_INPUT) ;
    }
    if (n < 0)      // or (n <= 0) ?
    {
        if (print_level > 0) {printf ("n invalid\n") ;}
        return (SLIP_INCORRECT_INPUT) ;
    }
    if (nzmax < 0)      // or (nzmax <= 0) ?
    {
        if (print_level > 0) {printf ("nzmax invalid\n") ;}
        return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------
    // Now check column/row indices
    //--------------------------------------------------------------------------
    switch (A->kind)
    {
        case SLIP_CSC:
        {
            int64_t* Ap = A->p;
            int64_t* Ai = A->i;

            //------------------------------------------------------------------
            // check the column pointers
            //------------------------------------------------------------------

            if (nz > 0 && (Ap == NULL || Ap [0] != 0 || Ap [n] != nz))
            {
                // column pointers invalid
                if (print_level > 0) {printf ("p invalid\n") ;}
                return (SLIP_INCORRECT_INPUT) ;
            }
            for (j = 0 ; j < n ; j++)
            {
                p = Ap [j] ;
                pend = Ap [j+1] ;
                if (pend < p || pend > nz)
                {
                    // column pointers not monotonically non-decreasing
                    if (print_level > 0) {printf ("p invalid\n") ;}
                    return (SLIP_INCORRECT_INPUT) ;
                }
            }

            //------------------------------------------------------------------
            // check the row indices && print values
            //------------------------------------------------------------------

            if (nzmax > 0 && (Ai == NULL || SLIP_X(A) == NULL))
            {
                // row indices or values not present
                if (print_level > 0) {printf ("i or x invalid\n") ;}
                return (SLIP_INCORRECT_INPUT) ;
            }

            // allocate workspace to check for duplicates
            mark = SLIP_calloc (n, sizeof (int64_t)) ;
            if (mark == NULL)
            {
                // out of memory
                if (print_level > 0) printf ("out of memory\n") ;
                SLIP_FREE_ALL;
                return (SLIP_OUT_OF_MEMORY) ;
            }

            for (j = 0 ; j < n ; j++)  // iterate across columns
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
                        if (print_level > 0) {printf ("invalid index\n") ;}
                        SLIP_FREE_ALL ;
                        return (SLIP_INCORRECT_INPUT) ;
                    }
                    if (print_level > 1)
                    {
                        printf ("  row %"PRId64" : ", i) ;
                        SLIP_info status = 0;

                        switch ( A->type)
                        {
                            case SLIP_MPZ:
                            {
                                status = SLIP_gmp_printf("%Zd \n", A->x.mpz[p]);
                                break;
                            }
                            case SLIP_MPQ:
                            {
                                status = SLIP_gmp_printf("%Qd \n", A->x.mpq[p]);
                                break;
                            }
                            case SLIP_MPFR:
                            {
                                status = SLIP_mpfr_printf("%.*Rf \n",
                                                          SLIP_GET_PRECISION(option),
                                                          A->x.mpfr [p]);
                                break;
                            }
                            case SLIP_FP64:
                            {
                                printf("%lf \n", A->x.fp64[p]);
                                break;
                            }
                            case SLIP_INT64:
                            {
                                printf("%ld \n", A->x.int64[p]);
                                break;
                            }
                            default: return SLIP_INCORRECT_INPUT;
                        }
                        if (status < 0)
                        {
                            SLIP_FREE_ALL ;
                            printf (" error: %d\n", status) ;
                            return (status) ;
                        }
                    }
                    mark [i] = marked ;
                }
            }
            break;
        }

        // Check triplet matrices
        case SLIP_TRIPLET:
        {

            int64_t* Aj = A->j;
            int64_t* Ai = A->i;

            //------------------------------------------------------------------
            // basic pointer checking
            //------------------------------------------------------------------

            if (nzmax > 0 && (Ai == NULL || Aj == NULL || SLIP_X(A) == NULL))
            {
                // row indices or values not present
                if (print_level > 0) printf ("i or j or x invalid\n") ;
                return (SLIP_INCORRECT_INPUT) ;
            }

            // allocate workspace to check for duplicates
            bmark = SLIP_calloc (m*n, sizeof (bool)) ;
            if (bmark == NULL)
            {
                // out of memory
                if (print_level > 0) printf ("out of memory\n") ;
                SLIP_FREE_ALL;
                return (SLIP_OUT_OF_MEMORY) ;
            }

            //------------------------------------------------------------------
            // check for duplicate and print each entry as "Ai Aj Ax"
            //------------------------------------------------------------------
            for (p = 0 ; p < nz ; p++)
            {
                i = Ai[p];
                j = Aj[p];
                if (i < 0 || i >= m || j < 0 || j >= n || bmark [j*m+i])
                {
                    // row indices out of range, or duplicate
                    if (print_level > 0) printf ("invalid index\n") ;
                    SLIP_FREE_ALL ;
                    return (SLIP_INCORRECT_INPUT) ;
                }
                if (print_level >= 2)
                {
                    printf ("  %"PRId64" %"PRId64" : ", i, j) ;
                    SLIP_info status = 0;

                    switch ( A->type)
                    {
                        case SLIP_MPZ:
                        {
                            status = SLIP_gmp_printf("%Zd \n", A->x.mpz [p]);
                            break;
                        }
                        case SLIP_MPQ:
                        {
                            status = SLIP_gmp_printf ("%Qd \n", A->x.mpq [p]);
                            break;
                        }
                        case SLIP_MPFR:
                        {
                            status = SLIP_mpfr_printf("%.*Rf \n", SLIP_GET_PRECISION(option),
                                                      A->x.mpfr [p]);
                            break;
                        }
                        case SLIP_FP64:
                        {
                            printf("%lf \n", A->x.fp64[p]);
                            break;
                        }
                        case SLIP_INT64:
                        {
                            printf("%ld \n", A->x.int64[p]);
                            break;
                        }
                        default: return SLIP_INCORRECT_INPUT;
                    }
                    if (status < 0)
                    {
                        SLIP_FREE_ALL ;
                        printf (" error: %d\n", status) ;
                        return (status) ;
                    }
                }
                bmark [i+j*m] = true ;
            }
            break;
        }
        // Check dense matrices
        case SLIP_DENSE:
        {
            // If A is dense, A->i, A->j etc are all NULL. All we must do is
            // to check that its dimensions are correct and print the values if
            // desired.

            if (nzmax > 0 && SLIP_X(A) == NULL)
            {
                // row indices or values not present
                if (print_level > 0) {printf ("x invalid\n") ;}
                return (SLIP_INCORRECT_INPUT) ;
            }
            //------------------------------------------------------------------
            // print values
            //------------------------------------------------------------------

            for (j = 0 ; j < n ; j++)
            {
                if (print_level >= 2)
                {
                    printf ("column %"PRId64" :\n", j) ;
                }
                for (i = 0; i < m; i++)
                {
                    if (print_level >= 2)
                    {
                        printf ("  row %"PRId64" : ", i) ;
                        SLIP_info status = 0;

                        switch ( A->type)
                        {
                            case SLIP_MPZ:
                            {
                                status = SLIP_gmp_printf ( "%Zd \n" ,
                                    SLIP_2D(A, i, j, mpz)) ;
                                break;
                            }
                            case SLIP_MPQ:
                            {
                                status = SLIP_gmp_printf ( "%Qd \n",
                                    SLIP_2D(A, i, j, mpq));
                                break;
                            }
                            case SLIP_MPFR:
                            {
                                status = SLIP_mpfr_printf ( "%.*Rf \n",
                                    SLIP_GET_PRECISION(option), SLIP_2D(A, i, j, mpfr));
                                break;
                            }
                            case SLIP_FP64:
                            {
                                printf("%lf \n", SLIP_2D(A, i, j, fp64));
                                break;
                            }
                            case SLIP_INT64:
                            {
                                printf("%ld \n", SLIP_2D(A, i, j, int64));
                                break;
                            }
                            default: return SLIP_INCORRECT_INPUT;
                        }
                        if (status < 0)
                        {
                            printf (" error: %d\n", status) ;
                            return (status) ;
                        }
                    }
                }
            }
            break;
        }
        default: return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE_ALL ;
    return (SLIP_OK) ;
}

