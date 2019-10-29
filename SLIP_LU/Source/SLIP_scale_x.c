# include "SLIP_LU_internal.h"

/* Purpose: This function scales the x matrix if necessary */

#define SLIP_FREE_WORKSPACE

SLIP_info SLIP_scale_x
(
    mpq_t **x,              // Solution matrix
    SLIP_sparse *A,         // matrix A
    SLIP_dense *b           // right hand side
)
{
    if (!x || !A || !b || !(A->scale) || !(b->scale)) 
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t r, n, numRHS;
    SLIP_info ok;
    n = A->m;
    numRHS = b->n;

    SLIP_CHECK(slip_mpq_cmp_ui(&r, A->scale, 1, 1));
    if (r != 0)
    {
        for (int32_t i = 0; i < n; i++)
        {
            for (int32_t j = 0; j < numRHS; j++)
            {
                SLIP_CHECK(slip_mpq_mul(x[i][j], x[i][j], A->scale));
            }
        }
    }

    SLIP_CHECK(slip_mpq_cmp_ui(&r, b->scale, 1, 1));
    if (r != 0)
    {
        for (int32_t i = 0; i < n; i++)
        {
            for (int32_t j = 0; j < numRHS; j++)
            {
                SLIP_CHECK(slip_mpq_div(x[i][j], x[i][j], b->scale));
            }
        }
    }
    return SLIP_OK;
}
#undef SLIP_FREE_WORKSPACE
