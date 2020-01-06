//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_dense: create an empty dense matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: Create an empty SLIP_dense matrix of size 0.
 * If there's no memory, instead returns NULL */

#define CHECK_RESULT(method)    \
{                               \
    if (method != SLIP_OK)      \
    {                           \
        SLIP_delete_dense(&A);  \
        return NULL;            \
    }                           \
}

SLIP_dense *SLIP_create_dense( void )
{
    SLIP_dense *A = SLIP_malloc(sizeof(SLIP_dense));
    // Check for out of memory
    if (A == NULL) { return NULL; }
    
    // Set m, n, and scale
    A->m = 0;
    A->n = 0;
    A->x = NULL;
    SLIP_MPQ_SET_NULL(A->scale);
    CHECK_RESULT (SLIP_mpq_init(A->scale));
    
    // Initial scale is 1
    CHECK_RESULT (SLIP_mpq_set_ui(A->scale, 1, 1));
    
    return A;
}

