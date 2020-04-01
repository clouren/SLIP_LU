//------------------------------------------------------------------------------
// SLIP_LU/SLIP_matrix_nnz: find # of entries in a matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "SLIP_LU_internal.h"

int64_t SLIP_matrix_nnz     // return # of entries in A, or -1 on error
(
    SLIP_matrix *A,         // matrix to query
    SLIP_options *option    // command options, currently unused
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (A == NULL ) 
    {
        return (-1) ;
    }

    //--------------------------------------------------------------------------
    // find nnz (A)
    //--------------------------------------------------------------------------

    switch (A->kind)
    {
        case SLIP_CSC:     return (A->p == NULL ? (-1) : A->p [A->n]) ;
        case SLIP_TRIPLET: return (A->nz) ;
        case SLIP_DENSE:   return (A->m * A->n) ;
        default:           return (-1) ;
    }
}

