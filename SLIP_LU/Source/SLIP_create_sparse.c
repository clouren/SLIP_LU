//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_sparse: create an empty sparse mpq matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function return a created empty sparse matrix as
 * SLIP_sparse pointer upon successful malloc
 */

#define CHECK_RESULT(method)    \
{                               \
    if (method != SLIP_OK)      \
    {                           \
        SLIP_delete_sparse(&A); \
        return NULL;            \
    }                           \
}

SLIP_sparse *SLIP_create_sparse( void )
{
    SLIP_sparse *A = (SLIP_sparse*) SLIP_malloc(sizeof(SLIP_sparse));
    if (A == NULL) {return NULL;}
    A -> n     = 0;
    A -> m     = 0;
    A -> nzmax = 0;
    A -> nz    = 0;
    A -> p     = NULL;
    A -> i     = NULL;
    A -> x     = NULL;

    SLIP_MPQ_SET_NULL(A->scale);
    CHECK_RESULT (SLIP_mpq_init(A->scale));

    // Initial scale is 1
    CHECK_RESULT (SLIP_mpq_set_ui(A->scale, 1, 1));

    return A;
}
