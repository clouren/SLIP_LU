//------------------------------------------------------------------------------
// SLIP_LU/SLIP_matrix_copy: create a copy of a matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// SLIP_matrix_copy creates a SLIP_matrix C that is a modified copy of a
// SLIP_matrix A.  The new matrix C can have a different kind and type
// than A.

// The input matrix A is assumed to be valid.  It can be checked first with
// SLIP_matrix_check, if desired.  If the input matrix A is not valid, results
// are undefined.

#define SLIP_FREE_WORK                  \
    SLIP_matrix_free (&T, option) ;     \
    SLIP_matrix_free (&Y, option) ;     \
    SLIP_FREE (W) ;

#define SLIP_FREE_ALL                   \
    SLIP_FREE_WORK ;                    \
    SLIP_matrix_free (&C, option) ;

#include "slip_internal.h"

SLIP_info SLIP_matrix_copy
(
    SLIP_matrix **C_handle, // matrix to create (never shallow)
    // inputs, not modified:
    SLIP_kind C_kind,       // C->kind: CSC, triplet, or dense
    SLIP_type C_type,       // C->type: mpz_t, mpq_t, mpfr_t, int64_t, or double
    SLIP_matrix *A,         // matrix to make a copy of (may be shallow)
    const SLIP_options *option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    int64_t nz = SLIP_matrix_nnz (A, option) ;
    if (A == NULL || C_handle == NULL || nz < 0 ||
        A->kind < SLIP_CSC || A->kind > SLIP_DENSE ||
        A->type < SLIP_MPZ || A->type > SLIP_FP64  ||
        C_kind  < SLIP_CSC || C_kind  > SLIP_DENSE ||
        C_type  < SLIP_MPZ || C_type  > SLIP_FP64)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }
    (*C_handle) = NULL ;
    SLIP_matrix *C = NULL ;
    SLIP_matrix *Y = NULL ;
    SLIP_matrix *T = NULL ;
    int64_t *W = NULL ;
    int64_t m = A->m ;
    int64_t n = A->n ;
    mpfr_rnd_t round = SLIP_OPTION_ROUND (option) ;

    //--------------------------------------------------------------------------
    // copy and convert A into C
    //--------------------------------------------------------------------------

    switch (C_kind)
    {

        //----------------------------------------------------------------------
        // C is CSC
        //----------------------------------------------------------------------

        case SLIP_CSC:
        {

            switch (A->kind)
            {

                //--------------------------------------------------------------
                // A is CSC, C is CSC
                //--------------------------------------------------------------

                case SLIP_CSC:
                {
                    // allocate C
                    SLIP_CHECK (SLIP_matrix_allocate (&C, SLIP_CSC, C_type,
                        m, n, nz, false, true, option)) ;
                    // copy the pattern of A into C
                    memcpy (C->p, A->p, (n+1) * sizeof (int64_t)) ;
                    memcpy (C->i, A->i, nz * sizeof (int64_t)) ;
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (SLIP_X (C), C->type,
                        SLIP_X (A), A->type, nz, C->scale, A->scale, option)) ;
                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is CSC
                //--------------------------------------------------------------

                case SLIP_TRIPLET:
                {

                    // Y = typecast the values of A into the type of C
                    // (not the pattern; Y is SLIP_DENSE)
                    SLIP_CHECK (slip_cast_matrix (&Y, C_type, A, option)) ;

                    // allocate workspace
                    W = (int64_t *) SLIP_calloc (n, sizeof (int64_t)) ;
                    if (W == NULL)
                    {
                        SLIP_FREE_ALL ;
                        return (SLIP_OUT_OF_MEMORY) ;
                    }

                    // allocate C
                    SLIP_CHECK (SLIP_matrix_allocate (&C, SLIP_CSC,
                        C_type, m, n, nz, false, true, option)) ;

                    // count the # of entries in each column
                    for (int64_t k = 0 ; k < nz ; k++)
                    {
                        W [A->j [k]]++ ;
                    }

                    // C->p = cumulative sum of W
                    slip_cumsum (C->p, W, n) ;

                    // build the matrix
                    switch (C->type)
                    {
                        case SLIP_MPZ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SLIP_CHECK (SLIP_mpz_set (
                                    SLIP_1D (C, p, mpz),
                                    SLIP_1D (Y, k, mpz))) ;
                            }
                            break ;

                        case SLIP_MPQ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SLIP_CHECK (SLIP_mpq_set (
                                    SLIP_1D (C, p, mpq),
                                    SLIP_1D (Y, k, mpq))) ;
                            }
                            break ;

                        case SLIP_MPFR:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SLIP_CHECK (SLIP_mpfr_set (
                                    SLIP_1D (C, p, mpfr),
                                    SLIP_1D (Y, k, mpfr),
                                    round)) ;
                            }
                            break ;

                        case SLIP_INT64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SLIP_1D (C, p, int64) =
                                    SLIP_1D (Y, k, int64) ;
                            }
                            break ;

                        case SLIP_FP64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SLIP_1D (C, p, fp64) =
                                    SLIP_1D (Y, k, fp64) ;
                            }
                            break ;

                    }

                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is CSC
                //--------------------------------------------------------------

                case SLIP_DENSE:
                {
                    // Y = typecast the values of A into the type of C
                    SLIP_CHECK (slip_cast_matrix (&Y, C_type, A, option)) ;
                    int s ;

                    // count the actual nonzeros in Y
                    if (Y->type == SLIP_FP64)
                    {
                    printf("\nhey nz is: %ld", nz);
                    printf("\nY->n and Y->m are: %ld %ld", Y->n, Y->m);
                    }
                    int64_t actual = 0 ;
                    switch (Y->type)
                    {

                        case SLIP_MPZ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                SLIP_CHECK (SLIP_mpz_sgn (&s,
                                    SLIP_1D (Y, k, mpz))) ;
                                if (s != 0) actual++ ;
                            }
                            break ;

                        case SLIP_MPQ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                SLIP_CHECK (SLIP_mpq_sgn (&s,
                                    SLIP_1D (Y, k, mpq))) ;
                                if (s != 0) actual++ ;
                            }
                            break ;

                        case SLIP_MPFR:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                SLIP_CHECK (SLIP_mpfr_sgn (&s,
                                    SLIP_1D (Y, k, mpfr))) ;
                                if (s != 0) actual++ ;
                            }
                            break ;

                        case SLIP_INT64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                if (SLIP_1D (Y, k, int64) != 0) actual++ ;
                            }
                            break ;

                        case SLIP_FP64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                if (SLIP_1D (Y, k, fp64) != 0) actual++ ;
                            }
                            break ;

                    }
                    if (Y->type == SLIP_FP64)
                    {
                        printf("\nHey actual is: %ld",actual);
                        printf("\nY[1,2] is %f", SLIP_2D(Y, 1, 2, fp64));
                        printf("\nThe index accessed is %ld", 1 + 2*Y->m);
                        printf("\nThe value in SLIP_1D is %f", Y->x.fp64[1 + 2*Y->m]);
                    }
                    // allocate C
                    SLIP_CHECK (SLIP_matrix_allocate (&C, SLIP_CSC, C_type,
                        m, n, actual, false, true, option)) ;

                    // Construct C
                    nz = 0 ;
                    switch (C->type)
                    {

                        case SLIP_MPZ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    //SLIP_CHECK (SLIP_mpz_sgn (&s,
                                    //    SLIP_2D (Y, i, j, mpz))) ;
                                    SLIP_CHECK( SLIP_mpz_sgn( &s, Y->x.mpz[ i + j*A->m]));
                                    if (s != 0)
                                    {
                                        C->i [nz] = i ;
                                        //SLIP_CHECK (SLIP_mpz_set (
                                        //    SLIP_1D (C, nz, mpz),
                                        //    SLIP_2D (Y, i, j, mpz))) ;
                                        SLIP_CHECK( SLIP_mpz_set ( SLIP_1D (C, nz, mpz),
                                                                   Y->x.mpz[ i + j*A->m] ));
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SLIP_MPQ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    //SLIP_CHECK (SLIP_mpq_sgn (&s,
                                    //    SLIP_2D (Y, i, j, mpq))) ;
                                    SLIP_CHECK (SLIP_mpq_sgn (&s,
                                        Y->x.mpq[ i + j*A->m])) ;
                                    if (s != 0)
                                    {
                                        C->i [nz] = i ;
                                        //SLIP_CHECK (SLIP_mpq_set (
                                        //    SLIP_1D (C, nz, mpq),
                                        //    SLIP_2D (Y, i, j, mpq))) ;
                                        SLIP_CHECK(SLIP_mpq_set (
                                            SLIP_1D(C, nz, mpq),
                                            Y->x.mpq[ i + j*A->m]));
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SLIP_MPFR:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    //SLIP_CHECK (SLIP_mpfr_sgn (&s,
                                    //    SLIP_2D (Y, i, j, mpfr))) ;
                                    SLIP_CHECK (SLIP_mpfr_sgn (&s,
                                        Y->x.mpfr[i + j*A->m])) ;
                                    if (s != 0)
                                    {
                                        C->i [nz] = i ;
                                        //SLIP_CHECK (SLIP_mpfr_set (
                                        //    SLIP_1D (C, nz, mpfr),
                                        //    SLIP_2D (Y, i, j, mpfr),
                                        //    round)) ;
                                        
                                        SLIP_CHECK (SLIP_mpfr_set (
                                            SLIP_1D (C, nz, mpfr),
                                            Y->x.mpfr[i + j*A->m],
                                            round)) ;
                                        
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SLIP_INT64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    //if (SLIP_2D (Y, i, j, int64) != 0)
                                    if ( Y->x.int64[i +j*A->m] != 0)
                                    {
                                        C->i [nz] = i ;
                                        //SLIP_1D (C, nz, int64) =
                                        //    SLIP_2D (Y, i, j, int64) ;
                                        SLIP_1D (C, nz, int64) =
                                            Y->x.int64[i +j*A->m] ;
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SLIP_FP64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    //if (SLIP_2D (Y, i, j, fp64) != 0)
                                    if ( Y->x.fp64[i +j*A->m] != 0)
                                    {
                                        C->i [nz] = i ;
                                        //SLIP_1D (C, nz, fp64) =
                                        //    SLIP_2D (Y, i, j, fp64) ;
                                        SLIP_1D (C, nz, fp64) =
                                            Y->x.fp64[i +j*A->m];
                                        nz++ ;
                                    }
                                    else
                                    {
                                        printf("\nY [%ld %ld] is zero", i,j);
                                        printf("\nReal value is %f", SLIP_2D(Y, i, j, fp64));
                                    }
                                }
                            }
                            break ;
                    }
                    C->p [n] = nz ;
                    printf("\nHey C->p[n] is: %ld", C->p[n]);
                }
                break ;

            }

        }
        break ;

        //----------------------------------------------------------------------
        // C is triplet
        //----------------------------------------------------------------------

        case SLIP_TRIPLET:
        {

            switch (A->kind)
            {

                //--------------------------------------------------------------
                // A is CSC, C is triplet
                //--------------------------------------------------------------

                case SLIP_CSC:
                {
                    // allocate C
                    SLIP_CHECK (SLIP_matrix_allocate (&C, SLIP_TRIPLET, C_type,
                        m, n, nz, false, true, option)) ;
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (SLIP_X (C), C->type,
                        SLIP_X (A), A->type, nz, C->scale, A->scale, option)) ;
                    // copy the row indices A->i into C->i
                    memcpy (C->i, A->i, nz * sizeof (int64_t)) ;
                    // construct C->j
                    for (int64_t j = 0 ; j < n ; j++)
                    {
                        for (int64_t p = A->p [j] ; p < A->p [j+1] ; p++)
                        {
                            C->j [p] = j ;
                        }
                    }
                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is triplet
                //--------------------------------------------------------------

                case SLIP_TRIPLET:
                {
                    // allocate C
                    SLIP_CHECK (SLIP_matrix_allocate (&C, SLIP_TRIPLET, C_type,
                        m, n, nz, false, true, option)) ;
                    // copy the pattern of A into C
                    memcpy (C->j, A->j, nz * sizeof (int64_t)) ;
                    memcpy (C->i, A->i, nz * sizeof (int64_t)) ;
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (SLIP_X (C), C->type,
                        SLIP_X (A), A->type, nz, C->scale, A->scale, option)) ;
                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is triplet
                //--------------------------------------------------------------

                case SLIP_DENSE:
                {
                    // convert A to a temporary CSC matrix
                    SLIP_CHECK (SLIP_matrix_copy (&T, SLIP_CSC, C_type,
                        A, option)) ;
                    // convert T from CSC to triplet
                    SLIP_CHECK (SLIP_matrix_copy (&C, SLIP_TRIPLET, C_type,
                        T, option)) ;
                    SLIP_matrix_free (&T, option) ;
                }
                break ;

            }

        }
        break ;

        //----------------------------------------------------------------------
        // C is dense
        //----------------------------------------------------------------------

        case SLIP_DENSE:
        {

            // allocate C
            SLIP_CHECK (SLIP_matrix_allocate (&C, SLIP_DENSE, C_type,
                m, n, nz, false, true, option)) ;

            switch (A->kind)
            {

                //--------------------------------------------------------------
                // A is CSC, C is dense
                //--------------------------------------------------------------

                case SLIP_CSC:
                {
                    // Y = typecast the values of A into the type of C
                    SLIP_CHECK (slip_cast_matrix (&Y, C->type, A, option)) ;

                    switch (C->type)
                    {

                        case SLIP_MPZ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SLIP_CHECK (SLIP_mpz_set (
                                        SLIP_2D (C, i, j, mpz),
                                        SLIP_1D (Y, p, mpz))) ;
                                }
                            }
                            break ;

                        case SLIP_MPQ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SLIP_CHECK (SLIP_mpq_set (
                                        SLIP_2D (C, i, j, mpq),
                                        SLIP_1D (Y, p, mpq))) ;
                                }
                            }
                            break ;

                        case SLIP_MPFR:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SLIP_CHECK (SLIP_mpfr_set (
                                        SLIP_2D (C, i, j, mpfr),
                                        SLIP_1D (Y, p, mpfr),
                                        round)) ;
                                }
                            }
                            break ;

                        case SLIP_INT64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SLIP_2D (C, i, j, int64) =
                                        SLIP_1D (Y, p, int64) ;
                                }
                            }
                            break ;

                        case SLIP_FP64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SLIP_2D (C, i, j, fp64) =
                                        SLIP_1D (Y, p, fp64) ;
                                }
                            }
                            break ;

                    }

                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is dense
                //--------------------------------------------------------------

                case SLIP_TRIPLET:
                {
                    // Y = typecast the values of A into the type of C
                    SLIP_CHECK (slip_cast_matrix (&Y, C->type, A, option)) ;

                    switch (C->type)
                    {

                        case SLIP_MPZ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SLIP_CHECK (SLIP_mpz_set (
                                    SLIP_2D (C, i, j, mpz),
                                    SLIP_1D (Y, k, mpz))) ;
                            }
                            break ;

                        case SLIP_MPQ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SLIP_CHECK (SLIP_mpq_set (
                                    SLIP_2D (C, i, j, mpq),
                                    SLIP_1D (Y, k, mpq))) ;
                            }
                            break ;

                        case SLIP_MPFR:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SLIP_CHECK (SLIP_mpfr_set (
                                    SLIP_2D (C, i, j, mpfr),
                                    SLIP_1D (Y, k, mpfr),
                                    round)) ;
                            }
                            break ;

                        case SLIP_INT64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SLIP_2D (C, i, j, int64) =
                                    SLIP_1D (Y, k, int64) ;
                            }
                            break ;

                        case SLIP_FP64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SLIP_2D (C, i, j, fp64) =
                                    SLIP_1D (Y, k, fp64) ;
                            }
                            break ;

                    }
                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is dense
                //--------------------------------------------------------------

                case SLIP_DENSE:
                {
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (SLIP_X (C), C->type,
                        SLIP_X (A), A->type, nz, C->scale, A->scale, option)) ;
                }
                break ;

            }

        }
        break ;

    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE_WORK ;
    (*C_handle) = C ;

    return (SLIP_OK) ;
}

