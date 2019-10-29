//------------------------------------------------------------------------------
// SLIP_LU/slip_sort_xi: sort the xi vector using the current row permutation
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function sorts the xi vector with respect to the current row
 * permutation. This sort is efficient as its complexity is |x| log |x|.
 * The idea of the sort is that you have xi[top, top+1, ...]. We essentially
 * mask them and sort the masked vector (which is with respect to the row
 * permutation). We then unmask them to get the correct value. For instance, the
 * correct sorted order could be [4 2 5 1] because of the column permutation.
 */

static inline int32_t compare (const void * a, const void * b)
{
    return ( *(int32_t*)a - *(int32_t*)b );
}

void slip_sort_xi
(
    int32_t* xi,        // nonzero pattern
    int32_t top,        // nonzeros are stored in xi[top..n-1]
    int32_t n,          // size of problem
    int32_t* pinv,      // inverse row permutation
    int32_t* row_perm   // opposite of pinv. if pinv[j] = k then row_perm[k] = j
)
{
    // Convert xi vector with respect to pinv
    for (int32_t j = top; j < n; j++)    
    {
        xi[j] = pinv[xi[j]];
    }
    // Sort xi[top..n-1]
    qsort(&xi[top], n-top, sizeof(int32_t), compare); 
    // Place xi back in original value
    for (int32_t j = top; j < n; j++)    
    {
        xi[j] = row_perm[xi[j]];
    }
}
