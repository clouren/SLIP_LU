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

#include "slip_internal.h"

/* check the validity of a SLIP_matrix */

// print_level from option struct:
//      0: nothing
//      1: just errors
//      2: errors and terse output
//      3: verbose

// TODO allow redirection of printf

SLIP_info SLIP_matrix_check     // returns a SLIP_LU status code
(
    const SLIP_matrix *A,     // matrix to check
    const SLIP_options* option
)
{

    //--------------------------------------------------------------------------
    // check the dimensions
    //--------------------------------------------------------------------------

    int print_level = SLIP_OPTION_PRINT_LEVEL(option);
    int64_t nz = SLIP_matrix_nnz(A, option);    // Number of nonzeros in A
    if (nz < 0)
    {
        if (print_level > 0) {printf ("A is NULL or nz invalid\n") ;}
        return (SLIP_INCORRECT_INPUT) ;
    }

    int64_t m = A->m ;
    int64_t n = A->n ;
    int64_t nzmax = A->nzmax ;

    if (m < 0)
    {
        if (print_level > 0) {printf ("m invalid\n") ;}
        return (SLIP_INCORRECT_INPUT) ;
    }
    if (n < 0)
    {
        if (print_level > 0) {printf ("n invalid\n") ;}
        return (SLIP_INCORRECT_INPUT) ;
    }
    if (nzmax < 0)
    {
        if (print_level > 0) {printf ("nzmax invalid\n") ;}
        return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------
    // check the dimensions
    //--------------------------------------------------------------------------

    if (A->type < SLIP_MPZ || A->type > SLIP_FP64)
    //  A->kind < SLIP_CSC || A->kind > SLIP_DENSE // checked in SLIP_matrix_nnz
    {
        if (print_level > 0)
        {
            printf ("A has invalid type.\n") ;
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
    // initialize workspace
    //--------------------------------------------------------------------------

    int64_t i, j, p, pend ;
    int64_t* mark = NULL;   // used for checking duplicate index of CSC
    bool* bmark = NULL;     // used for checking duplicate index of triplet
    uint64_t prec = SLIP_OPTION_PREC (option);

    //--------------------------------------------------------------------------
    // check each kind of matrix: CSC, triplet, or dense
    //--------------------------------------------------------------------------

    switch (A->kind)
    {

        //----------------------------------------------------------------------
        // check a matrix in CSC format
        //----------------------------------------------------------------------

        case SLIP_CSC:
        {
            int64_t* Ap = A->p;
            int64_t* Ai = A->i;

            //------------------------------------------------------------------
            // check the column pointers
            //------------------------------------------------------------------

            if (nzmax > 0 && (Ap == NULL || Ap [0] != 0))
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
            mark = (int64_t *) SLIP_calloc (m, sizeof (int64_t)) ;
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
                                    prec, A->x.mpfr [p]);
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
        }
        break;

        //----------------------------------------------------------------------
        // check a matrix in triplet format
        //----------------------------------------------------------------------

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
            bmark = (bool *) SLIP_calloc (m*n, sizeof (bool)) ;
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
                            status = SLIP_mpfr_printf("%.*Rf \n",
                                prec, A->x.mpfr [p]);
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
        }
        break;

        //----------------------------------------------------------------------
        // check a matrix in dense format
        //----------------------------------------------------------------------

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
                                    prec, SLIP_2D(A, i, j, mpfr));
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
                        }
                        if (status < 0)
                        {
                            printf (" error: %d\n", status) ;
                            return (status) ;
                        }
                    }
                }
            }
        }
        break;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE_ALL ;
    return (SLIP_OK) ;
}

