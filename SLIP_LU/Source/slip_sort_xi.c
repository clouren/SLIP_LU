//------------------------------------------------------------------------------
// SLIP_LU/slip_sort_xi: sort the xi vector using the current row permutation
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function sorts the xi vector with respect to the current row
 * permutation. This sort is efficient as its complexity is |x| log |x|.  The
 * idea of the sort is that you have xi[top, top+1, ...]. We essentially mask
 * them and sort the masked vector (which is with respect to the row
 * permutation). We then unmask them to get the correct value. For instance,
 * the correct sorted order could be [4 2 5 1] because of the column
 * permutation.
 */

#include "SLIP_LU_internal.h"

// TODO: this is used only by slip_ref_triangular_solve.
// Remove this function and put the code in that function.

static inline int compare (const void * a, const void * b)
{
    int64_t delta = ( *(int64_t*)a - *(int64_t*)b ) ;
    if (delta < 0)
    {
        return (-1) ;
    }
    else if (delta > 0)
    {
        return (1) ;
    }
    else
    {
        return (0) ;
    }
}

void slip_sort_xi
(
    int64_t* xi,              // nonzero pattern
    int64_t top,              // nonzeros are stored in xi[top..n-1]
    int64_t n,                // size of problem
    const int64_t* pinv,      // inverse row permutation
    const int64_t* row_perm   // opposite of pinv. if pinv[j] = k,
                              //    then row_perm[k] = j
)
{

    // Convert xi vector with respect to pinv
    for (int64_t j = top; j < n; j++)
    {
        xi[j] = pinv[xi[j]];
    }

    // Sort xi[top..n-1]
    qsort(&xi[top], n-top, sizeof(int64_t), compare);

    // Place xi back in original value
    for (int64_t j = top; j < n; j++)
    {
        xi[j] = row_perm[xi[j]];
    }
}

