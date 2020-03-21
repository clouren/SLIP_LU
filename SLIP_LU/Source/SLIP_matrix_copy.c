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

// TODO: this is a sketch (it won't compile yet).  See the TODO's.
// probably need some helper functions, since some chunks of code are
// repeated.

// The input matrix A is assumed to be valid.  It can be checked first with
// SLIP_matrix_check, if desired.  If the input matrix A is not valid, results
// are undefined.

// TODO: SLIP_matrix_check replaces SLIP_spok and works for any kind and type.

//  What to do about C->scale ?
//      A (double) to C (in mpz).  C->scale = ... ?
//      A (mpz) to C (in double).  use A->scale ?
//      I guess this function needs to compute C->scale somehow.
//      Does the user need to specify anything?

#define SLIP_FREE_WORK                  \
    free Y_mpq /* TODO */ \
    free Y_mpz /* TODO */ \
    free Y_mpfr /* TODO */ \
    SLIP_FREE (Y_int32) ;               \
    SLIP_FREE (Y_fp64) ;                \
    SLIP_FREE (W) ;

#define SLIP_FREE_ALL                   \
    SLIP_FREE_WORK ;                    \
    SLIP_matrix_free (&C, option) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_matrix_copy
(
    SLIP_matrix **C_handle, // matrix to create (never shallow)
    // inputs, not modified:
    SLIP_kind kind,         // CSC, triplet, or dense
    SLIP_type type,         // mpz_t, mpq_t, mprf_t, int32_t, or double
    SLIP_matrix *A,         // matrix to make a copy of (may be shallow)
    SLIP_options *option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    if (A == NULL || option == NULL || C_handle == NULL)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }
    (*C_handle) = NULL ;
    SLIP_matrix *C = NULL ;
    int32_t *W = NULL ;
    int32_t m = A->m ;
    int32_t n = A->n ;

    // workspace for each data type
    mpq_t   *Y_mpq   = NULL ;
    mpz_t   *Y_mpz   = NULL ;
    mpfr_t  *Y_mpfr  = NULL ;
    int32_t *Y_int32 = NULL ;
    double  *Y_fp64  = NULL ;

#if THIS_IS_A_DRAFT_AND_WONT_COMPILE_YET

    //--------------------------------------------------------------------------
    // copy and convert A into C
    //--------------------------------------------------------------------------

    switch (kind)
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
                    int32_t nz = A->p [n] ;
                    SLIP_CHECK (SLIP_matrix_allocate (&C, kind, type, m, n, nz,
                        false, option)) ;
                    // copy the pattern of A into C
                    memcpy (C->p, A->p, (n+1) * sizeof (int32)) ;
                    memcpy (C->i, A->i, nz * sizeof (int32)) ;
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (C->x, type, A->x, A->type,
                        nz, C->scale, option)) ;
                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is CSC
                //--------------------------------------------------------------

                case SLIP_TRIPLET:
                {
                    int32_t nz = A->nz ;

                    // Y = typecasted/scaled copy of A->x
                    allocate Y with the right type, of size nz  // TODO
                    SLIP_CHECK (slip_cast_array (Y, type, A->x, A->type,
                        nz, scale, option)) ;

                    // allocate workspace
                    W = (int32_t *) SLIP_calloc (n, sizeof (int32_t)) ;
                    if (W == NULL)
                    {
                        SLIP_FREE_ALL ;
                        return (SLIP_OUT_OF_MEMORY) ;
                    }
                    // allocate C
                    SLIP_CHECK (SLIP_matrix_allocate (&C, kind, type, m, n, nz,
                        false, option)) ;
                    C->scale = scale ; // TODO
                    // count the # of entries in each column
                    for (int32_t k = 0 ; k < nz ; k++)
                    {
                        W [A->j [k]]++ ;
                    }
                    // C->p = cumulative sum of W
                    slip_cumsum (C->p, W, n) ;
                    // build the matrix
                    // TODO need 5-way switch for typecast
                    // (Y and C already have the same type)
                    for (int32_t k = 0 ; k < nz ; k++)
                    {
                        // find the position in C for the kth tuple
                        int32_t p = W [J [k]]++ ;
                        // place the entry in C
                        C->i [p] = I [k] ;
                        SLIP_ENTRY (C, p, type) = SLIP_ENTRY (Y, k, type) ;
                    }
                    C->nz = C->p [n] ;          // TODO is this needed?
                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is CSC
                //--------------------------------------------------------------

                case SLIP_DENSE:
                {
                    int32_t nzmax = m*n ;
                    // Y = typecasted/scaled copy of A->x
                    allocate Y with the right type, of size nzmax   // TODO
                    SLIP_CHECK (slip_cast_array (Y, type, A->x, A->type,
                        nzmax, scale, option)) ;
                    // count the nonzeros in Y
                    int32_t nz = 0 ;
                    // TODO: make this a helper function with 5-way switch:
                    for (int32_t k = 0 ; k < nzmax ; k++)
                    {
                        if (SLIP_ENTRY (Y, k, type) != 0)  // TODO
                        {
                            nz++ ;
                        }
                    }
                    // allocate C
                    SLIP_CHECK (SLIP_matrix_allocate (&C, kind, type, m, n, nz,
                        false, option)) ;
                    C->scale = scale ; // TODO
                    // Construct C
                    // TODO: need 5-way switch for each type:
                    // (Y and C already have the same type)
                    int32_t nz = 0 ;
                    for (int32_t j = 0 ; j < n ; j++)
                    {
                        C->p [j] = nz ;
                        for (int32_t i = 0 ; i < m ; i++)
                        {
                            if (SLIP_2D (Y, i, j, type) != 0)  // TODO
                            {
                                C->i [nz] = i ;
                                SLIP_ENTRY (C, nz, type) =
                                    SLIP_2D (Y, i, j, type) ;  // TODO
                                nz++ ;
                            }
                        }
                    }
                    C->p [n] = nz ;
                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
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
                    int32_t nz = A->p [n] ;
                    SLIP_CHECK (SLIP_matrix_allocate (&C, kind, type, m, n, nz,
                        false, option)) ;
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (C->x, type, A->x, A->type,
                        nz, C->scale, option)) ;
                    // copy the row indices A->i into C->i
                    memcpy (C->i, A->i, nz * sizeof (int32)) ;
                    // construct C->j
                    for (int32_t j = 0 ; j < n ; j++)
                    {
                        for (int32_t p = A->p [j] ; p < A->p [j+1] ; p++)
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
                    int32_t nz = A->nz ;
                    SLIP_CHECK (SLIP_matrix_allocate (&C, kind, type, m, n, nz,
                        false, option)) ;
                    // copy the pattern of A into C
                    memcpy (C->j, A->j, nz * sizeof (int32)) ;
                    memcpy (C->i, A->i, nz * sizeof (int32)) ;
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (C->x, type, A->x, A->type,
                        nz, C->scale, option)) ;
                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is triplet
                //--------------------------------------------------------------

                case SLIP_DENSE:
                {
                    // convert A to a temporary CSC matrix, then to triplet
                    SLIP_CHECK (SLIP_matrix_copy (&T, SLIP_CSC, type, A,
                        option)) ;
                    SLIP_CHECK (SLIP_matrix_copy (&C, kind, type, T, option)) ;
                    SLIP_matrix_free (&T, option) ;

                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
            }

        }
        break ;

        //----------------------------------------------------------------------
        // C is dense
        //----------------------------------------------------------------------

        case SLIP_DENSE:
        {

            // allocate C
            int32_t nzmax = m*n ;
            SLIP_CHECK (SLIP_matrix_allocate (&C, kind, type, m, n, nzmax,
                false, option)) ;

            switch (A->kind)
            {

                //--------------------------------------------------------------
                // A is CSC, C is dense
                //--------------------------------------------------------------

                case SLIP_CSC:
                {
                    // TODO: need 5-way switch for each type:
                    for (int32_t j = 0 ; j < n ; k++)
                    {
                        for (int32_t p = A->p [j] ; p < A->p [j+1] ; p++)
                        {
                            int32_t i = A->i [p] ;
                            SLIP_2D (C, i, j, type) = SLIP_ENTRY (Y, p, type) ;
                        }
                    }
                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is dense
                //--------------------------------------------------------------

                case SLIP_TRIPLET:
                {
                    // TODO: need 5-way switch for each type:
                    int32_t nz = A->nz ;
                    for (int32_t k = 0 ; k < nz ; k++)
                    {
                        int32_t i = A->i [k] ;
                        int32_t j = A->j [k] ;
                        SLIP_2D (C, i, j, type) = A->x [k] ;    // TODO
                    }
                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is dense
                //--------------------------------------------------------------

                case SLIP_DENSE:
                {
                    // copy and typecast A->x into C->x
                    SLIP_CHECK (slip_cast_array (C->x, type, A->x, A->type,
                        nzmax, C->scale, option)) ;
                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
            }

        }
        break ;

        default: return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SLIP_FREE_WORK ;
    (*C_handle) = C ;
#endif

    return (SLIP_OK) ;
}

