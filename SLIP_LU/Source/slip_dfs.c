//------------------------------------------------------------------------------
// SLIP_LU/slip_dfs: depth-first search
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function performs a depth first search of the graph of the
 * matrix starting at node j. The output of this function is the set of nonzero
 * indices in the xi vector
 * 
 * This function is modified from CSparse/cs_dfs.
 */
void slip_dfs // performs a dfs of the graph of the matrix starting at node j
(
    int32_t *top,    // beginning of stack
    int32_t j,       // What node to start DFS at
    SLIP_sparse* L,  // matrix which represents the Graph of L 
    int32_t* xi,     // the nonzero pattern
    int32_t* pstack, // workspace vector
    int32_t* pinv    // row permutation 
)
{
    // No check here, input is checked in slip_reach.c
    int32_t i, p, p2, done, jnew, head = 0;
    
    // Initialize the recursion stack
    xi[0] = j;

    while (head >= 0)
    {
        // The j value of the nonzero
        j = xi[head];
        // The relative j value 
        jnew = pinv[j];

        //----------------------------------------------------------------------
        // Mark L->p[j] if not marked yet
        //----------------------------------------------------------------------
        if (!SLIP_MARKED (L->p,j))
        {
            SLIP_MARK(L->p,j);
            pstack[head] = (jnew < 0) ? 0 : SLIP_UNFLIP(L->p[jnew]);
        }
        // Node j is done if no unvisited neighbors
        done = 1;

        p2 = (jnew < 0) ? 0 : SLIP_UNFLIP(L->p[jnew+1]);

        //----------------------------------------------------------------------
        // Examine all neighbors of j
        //----------------------------------------------------------------------
        for (p = pstack[head]; p < p2; p++)
        {
            // Looking at neighbor node i
            i = L->i[p];
            // Skip already visited node
            if (SLIP_MARKED(L->p,i))  {continue;}
            
            // pause DFS of node j
            pstack[head] = p;
            // Start DFS at node i
            xi[++head] = i;
            // node j is not done
            done = 0;
            // break to start dfs i
            break;
        }
        if (done != 0)
        {
            head--;
            xi[--(*top)] = j;
        }
    }
    return;
}
