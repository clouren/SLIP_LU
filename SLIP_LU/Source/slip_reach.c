//------------------------------------------------------------------------------
// SLIP_LU/slip_reach: compute the set of nodes reachable from an input set
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

/* 
 * Purpose: This function computes the reach of column k of A on the graph of L
 * mathematically that is: xi = Reach(A(:,k))_G_L
 *
 * This function is derived from CSparse/cs_reach.c
 */
void slip_reach    // compute the reach of column k of A on the graph of L
(
    int32_t *top,
    SLIP_sparse* L,   // matrix representing graph of L
    SLIP_sparse* A,   // input matrix 
    int32_t k,        // column of A of interest
    int32_t* xi,      // nonzero pattern
    int32_t* pinv     // row permutation
)
{
    // inputs have been checked in slip_REF_triangular_solve
    int32_t p, n = L->n;
    *top = n;

    //--------------------------------------------------------------------------
    // Iterating across number of nonzero in column k
    //--------------------------------------------------------------------------
    for (p = A->p[k]; p < A->p[k + 1]; p++)
    { 
        // DFS at unmarked node i
        if (!SLIP_MARKED(L->p, A->i[p]))
        {
            slip_dfs(top, A->i[p], L, xi, xi+n, pinv);
        }
    }
    
    // Restore L
    for ( p = *top; p < n; p++)        
    {
        SLIP_MARK(L->p, xi[p]);
    }
    return;
}
