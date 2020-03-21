//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_double_mat: delete double dense matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function frees a dense double matrix.
 *
 * Input is a double*** mat and is destroyed on function completion.
 */

#include "SLIP_LU_internal.h"

// ignore warnings about unused parameters in this file
#pragma GCC diagnostic ignored "-Wunused-parameter"

void SLIP_delete_double_mat
(
    double*** A,   // dense matrix
    int64_t m,     // number of rows of A
    int64_t n      // number of columns of A (unused parameter)
)
{

    // TODO: use SLIP_matrix_free (&A, option) ;
    // Delete this function since *_mat will no longer be used.

    if (A == NULL || (*A) == NULL) {return;}

    for (int64_t i = 0; i < m; i++)
    {
        SLIP_FREE( (*A)[i]);
    }
    SLIP_FREE(*A);
}

