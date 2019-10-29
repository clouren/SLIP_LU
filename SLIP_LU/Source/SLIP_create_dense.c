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
    if (A == NULL) { return A; }
    
    // Set m, n, and scale
    A->m = 0;
    A->n = 0;
    A->x = NULL;
    SLIP_MPQ_SET_NULL(A->scale);
    CHECK_RESULT (slip_mpq_init(A->scale));
    
    // Initial scale is 1
    CHECK_RESULT (slip_mpq_set_ui(A->scale, 1, 1));
    
    return A;
}
#undef CHECK_RESULT
