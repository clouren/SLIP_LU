//------------------------------------------------------------------------------
// SLIP_LU/SLIP_check_solution: check solution to Ax=b
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#define SLIP_FREE_WORKSPACE                 \
    SLIP_MPQ_CLEAR(temp);                   \
    SLIP_delete_mpq_mat(&b2, n, numRHS);

# include "SLIP_LU_internal.h"

/* ========================================================================== */
/* ============= Check the solution of the linear system===================== */
/* ========= Performs a very quick rational arithmetic A*x=b ================ */
/* ========================================================================== */

SLIP_info SLIP_check_solution
(
    SLIP_sparse *A,          // input matrix
    mpq_t** x,               // solution vector
    SLIP_dense *b            // righ hand side
)
{
    if (!A || !x || !b || !b->x || !A->p || !A->i || !A->x) 
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare vars
    //--------------------------------------------------------------------------

    int32_t p, j, i, r;
    int32_t n = A->n, numRHS = b->n;
    mpz_t** bx = b->x;
    SLIP_info ok;
    mpq_t** b2 = NULL;
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(SLIP_mpq_init(temp));
    
    // b2 stores the solution of A*x
    b2 = SLIP_create_mpq_mat(n, numRHS);
    if (!b2) 
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // perform mpq_addmul in loops
    //--------------------------------------------------------------------------

    for (j = 0; j < numRHS; j++)
    {
        for (i = 0; i < n; i++)
        {
            for (p = A->p[i]; p < A->p[i + 1]; p++)
            {
                // temp = A[p][i]
                SLIP_CHECK(SLIP_mpq_set_z(temp, A->x[p]));

                // temp = temp*x[i]
                SLIP_CHECK(SLIP_mpq_mul(temp, temp, x[i][j]));

                // b2[p] = b2[p]-temp
                SLIP_CHECK(SLIP_mpq_add(b2[A->i[p]][j], b2[A->i[p]][j],temp));
            }
        }
    }

    //--------------------------------------------------------------------------
    // check if b==b2
    //--------------------------------------------------------------------------
    
    for (j = 0; j < numRHS; j++)
    {
        for (i = 0; i < n; i++)
        {
            // z = b[i] (correct b)
            SLIP_CHECK(SLIP_mpq_set_z(temp, bx[i][j]));

            // set check false if b!=b2
            SLIP_CHECK(SLIP_mpq_equal(&r, temp, b2[i][j]));
            if (r == 0)
            {
                SLIP_FREE_WORKSPACE;
                return SLIP_INCORRECT;
            }
        }
    }
    
    //--------------------------------------------------------------------------
    // Free memory 
    //--------------------------------------------------------------------------

    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

