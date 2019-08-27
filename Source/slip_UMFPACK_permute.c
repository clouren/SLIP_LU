# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function permutes the b vector & A matrix for the UMFPACK
 * ordering. This is necessary since we use UMFPACK to preorder A, thus we
 * must properly maintain b as well.
 */

#define SLIP_FREE_WORKSPACE             \
    SLIP_FREE(pinv);                    \
    SLIP_delete_mpz_mat(&b2, n, numRHS);

SLIP_info slip_UMFPACK_permute
(
    SLIP_dense *b,   // RHS vectors
    SLIP_sparse* A,  // input matrix
    int32_t* p_umf   // row permutation given by UMFPACK 
)
{
    if(!b || !b->x || !p_umf || !A || !A->p || !A->i || !A->x) 
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t k, i, index, n = A->n, numRHS = b->n;
    mpz_t** bx = b->x;
    SLIP_info ok = SLIP_OK;
    int32_t* pinv = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
    mpz_t** b2 = SLIP_create_mpz_mat(n, numRHS);
    if (!pinv || !b2) 
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    for (k = 0; k < n; k++)
    {
        index = p_umf[k];
        pinv[index] = k;
    }
    for (k = 0; k < A->nz; k++)        
    {
        A->i[k] = pinv[A->i[k]];
    }
    for (i = 0; i < numRHS; i++)
    {
        for (k = 0; k < n; k++)
        {
            SLIP_CHECK(slip_mpz_set(b2[pinv[k]][i], bx[k][i]));
        }
    }
    for (i = 0; i < numRHS; i++)
    {
        for (k = 0; k < n; k++)
        {
            SLIP_CHECK(slip_mpz_set(bx[k][i], b2[k][i]));
        }
    }
    SLIP_FREE_WORKSPACE;
    return ok;
}
#undef SLIP_FREE_WORKSPACE
