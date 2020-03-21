//------------------------------------------------------------------------------
// SLIP_LU/SLIP_matrix_free: free a SLIP_matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Free a SLIP_matrix.  Any shallow component is not freed.

#include "SLIP_LU_internal.h"

SLIP_info SLIP_matrix_free
(
    SLIP_matrix **A_handle, // matrix to free
    SLIP_options *option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (option == NULL)
    {
        // option is currently unused, but it's checked anyway
        return (SLIP_INCORRECT_INPUT) ;
    }
    if (A_handle == NULL || (*A_handle) == NULL)
    {
        // nothing to free (not an error)
        return (SLIP_OK) ;
    }
    SLIP_matrix *A = (*A_handle) ;

    //--------------------------------------------------------------------------
    // free any non-shallow components
    //--------------------------------------------------------------------------

    // free the integer pattern
    if (!(A->p_shallow)) SLIP_FREE (A->p) ;
    if (!(A->i_shallow)) SLIP_FREE (A->i) ;
    if (!(A->j_shallow)) SLIP_FREE (A->j) ;

    // free the values
    if (!(A->x_shallow))
    {
        switch (A->type)
        {
            case SLIP_MPZ:
                SLIP_delete_mpz_array (&(A->x.mpz), A->nzmax) ;
                break ;

            case SLIP_MPQ:
                SLIP_delete_mpq_array (&(A->x.mpq), A->nzmax) ;
                break ;

            case SLIP_MPFR:
                SLIP_delete_mpfr_array (&(A->x.mpfr), A->nzmax) ;
                break ;

            case SLIP_INT64:
                SLIP_FREE (A->x.int64) ;
                break ;

            case SLIP_FP64:
                SLIP_FREE (A->x.fp64) ;
                break ;

            default:
                // do nothing
                break ;
        }
    }

    // A->scale is never shallow
    SLIP_MPQ_CLEAR (A->scale) ;

    //--------------------------------------------------------------------------
    // free the header
    //--------------------------------------------------------------------------

    // the header is never shallow
    SLIP_FREE (A) ;
    (*A_handle) = NULL ;
    return (SLIP_OK) ;
}

