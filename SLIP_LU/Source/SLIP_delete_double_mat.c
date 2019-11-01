//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_double_mat: delete double dense matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: This function frees a dense double matrix.
 *
 * Input is a double*** mat and is destroyed on function completion.
 */
void SLIP_delete_double_mat
(
    double*** A,   // dense matrix
    int32_t m,     // number of rows of A
    int32_t n      // number of columns of A
)
{
    if (A == NULL || (*A) == NULL) {return;}
    for (int32_t i = 0; i < m; i++)
    {
        SLIP_FREE( (*A)[i]);
    }
    SLIP_FREE(*A);
}

