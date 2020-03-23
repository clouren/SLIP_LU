//------------------------------------------------------------------------------
// SLIP_LU/SLIP_matrix_allocate: allocate a SLIP_matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Allocate an m-by-n SLIP_matrix, in one of 15 data structures:
// (sparse CSC, sparse triplet, or dense) x
// (mpz, mpz, mfpr, int64, or double).

// The matrix may be created as 'shallow', in which case A->p, A->i, A->j, and
// A->x are all returned as NULL, and all A->*_shallow flags are returned as
// true.

#define SLIP_FREE_ALL \
    SLIP_matrix_free (&A, option) ;

#include "SLIP_LU_internal.h"

SLIP_info SLIP_matrix_allocate
(
    SLIP_matrix **A_handle, // matrix to allocate
    SLIP_kind kind,         // CSC, triplet, or dense
    SLIP_type type,         // mpz, mpq, mpfr, int64, or double
    int64_t m,              // # of rows
    int64_t n,              // # of columns
    int64_t nzmax,          // max # of entries for CSC or triplet
    bool shallow,           // if true, matrix is shallow.  A->p, A->i, A->j,
                            // A->x are all returned as NULL and must be set
                            // by the caller.  All A->*_shallow are returned
                            // as true.
    // TODO add this?:
    // bool init,
    SLIP_options *option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    if (A_handle == NULL)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }
    (*A_handle) = NULL ;
    if (m < 0 || n < 0 || option == NULL)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------
    // allocate the header
    //--------------------------------------------------------------------------

    SLIP_matrix *A = (SLIP_matrix *) SLIP_calloc (1, sizeof (SLIP_matrix)) ;
    if (A == NULL)
    {
        return (SLIP_OUT_OF_MEMORY) ;
    }

    if (kind == SLIP_DENSE)
    {
        nzmax = m*n ;
    }
    nzmax = SLIP_MAX (nzmax, 1) ;

    A->m = m ;
    A->n = n ;
    A->nzmax = nzmax ;
    A->nz = 0 ;
    A->kind = kind ;
    A->type = type ;

    // A->p, A->i, A->j, and A->x are currently NULL since A was calloc'd.
    A->p_shallow = false ;
    A->i_shallow = false ;
    A->j_shallow = false ;
    A->x_shallow = false ;

    // A->scale = 1
    SLIP_MPQ_SET_NULL (A->scale) ;
    SLIP_CHECK (SLIP_mpq_init (A->scale)) ;
    SLIP_CHECK (SLIP_mpq_set_ui (A->scale, 1, 1)) ;

    //--------------------------------------------------------------------------
    // allocate the p, i, j, and x components
    //--------------------------------------------------------------------------

    if (shallow)
    {

        // all components are shallow.  The caller can modify individual
        // components after A is created, as needed.
        A->p_shallow = true ;
        A->i_shallow = true ;
        A->j_shallow = true ;
        A->x_shallow = true ;

    }
    else
    {

        bool ok = true ;

        // allocate the integer pattern
        switch (kind)
        {
            case SLIP_CSC:
                A->p = (int64_t *) SLIP_calloc (n+1, sizeof (int64_t)) ;
                A->i = (int64_t *) SLIP_malloc (nzmax * sizeof (int64_t)) ;
                ok = (A->p != NULL && A->i != NULL) ;
                break ;

            case SLIP_TRIPLET:
                A->i = (int64_t *) SLIP_malloc (nzmax * sizeof (int64_t)) ;
                A->j = (int64_t *) SLIP_malloc (nzmax * sizeof (int64_t)) ;
                ok = (A->i != NULL && A->j != NULL) ;
                break ;

            case SLIP_DENSE:
                // nothing to do
                break ;

            default:
                SLIP_FREE_ALL ;
                return (SLIP_INCORRECT_INPUT) ;
        }

        // allocate the values
        switch (type)
        {
            case SLIP_MPZ:
                A->x.mpz = SLIP_create_mpz_array (nzmax) ;  // TODO
                ok = ok && (A->x.mpz != NULL) ;
                break ;

            case SLIP_MPQ:
                A->x.mpq = SLIP_create_mpq_array (nzmax) ;  // TODO
                ok = ok && (A->x.mpq != NULL) ;
                break ;

            case SLIP_MPFR:
                A->x.mpfr = SLIP_create_mpfr_array (nzmax, option) ;    // TODO
                ok = ok && (A->x.mpfr != NULL) ;
                break ;

            case SLIP_INT64:
                A->x.int64 = (int64_t *) SLIP_calloc (nzmax, sizeof (int64_t)) ;
                ok = ok && (A->x.int64 != NULL) ;
                break ;

            case SLIP_FP64:
                A->x.fp64 = (double *) SLIP_calloc (nzmax, sizeof (double)) ;
                ok = ok && (A->x.fp64 != NULL) ;
                break ;

            default:
                SLIP_FREE_ALL ;
                return (SLIP_INCORRECT_INPUT) ;
        }

        if (!ok)
        {
            SLIP_FREE_ALL ;
            return (SLIP_INCORRECT_INPUT) ;
        }
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    (*A_handle) = A ;
    return (SLIP_OK) ;
}

