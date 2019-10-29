//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_LU_analysis: TODO what does this do??
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function returns a pointer to an created SLIP_col type
 * with the length of S->q set as n (which needs to 1 + number of rows of input
 * matrix) upon successful malloc, otherwise, return NULL
 */

SLIP_LU_analysis *SLIP_create_LU_analysis
(
    int32_t n     //length of S->q
)
{
    SLIP_LU_analysis* S = NULL;
    S = (SLIP_LU_analysis*) SLIP_malloc(sizeof(SLIP_LU_analysis));
    if (S == NULL) {return S;}
    
    S->q = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
    if (S->q == NULL) 
    {
	SLIP_FREE(S);
	return NULL;
    }
    
    S->lnz = 0;
    S->unz = 0;
    return S;
}
