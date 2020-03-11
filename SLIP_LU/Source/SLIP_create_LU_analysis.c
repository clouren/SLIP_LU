//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_LU_analysis: Allocate memory for symbolic analysis
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function allocates memory for the symbolic analysis struct
 */

#include "SLIP_LU_internal.h"

SLIP_LU_analysis *SLIP_create_LU_analysis
(
    int32_t n     // number of columns in matrix to be analyzed
)
{
    // ALlocate memory for S
    SLIP_LU_analysis* S = NULL;
    S = (SLIP_LU_analysis*) SLIP_malloc(sizeof(SLIP_LU_analysis));
    if (S == NULL) {return S;}

    // Allocate memory for column permutation
    S->q = (int32_t*) SLIP_malloc((n+1) * sizeof(int32_t));
    if (S->q == NULL)
    {
        SLIP_FREE(S);
        return NULL;
    }

    S->lnz = 0;
    S->unz = 0;
    return S;
}
