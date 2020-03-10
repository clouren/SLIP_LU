//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_int_mat: delete a dense int matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function deletes a dense int matrix.
 *
 * Input is a int*** mat and its dimensions.
 * Input mat is destroyed on function completion
 */

#include "SLIP_LU_internal.h"

// ignore warnings about unused parameters in this file
#pragma GCC diagnostic ignored "-Wunused-parameter"

void SLIP_delete_int_mat
(
    int32_t ***A,  // dense matrix
    int32_t m,     // number of rows
    int32_t n      // number of columns (unused parameter)
)
{
    if (A == NULL || (*A) == NULL) {return;}
    for (int32_t i = 0; i < m; i++)
    {
        SLIP_FREE( (*A)[i]);
    }
    SLIP_FREE(*A);
}

