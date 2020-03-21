//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_mpq_mat: delete a dense mpq matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function deletes a dense mpq matrix.
 *
 * Input is a mpq_t*** matrix which is destroyed upon function completion.
 */

#include "SLIP_LU_internal.h"

void SLIP_delete_mpq_mat
(
    mpq_t***A,     // dense mpq matrix
    int64_t m,     // number of rows of A
    int64_t n      // number of columns of A
)
{

    // TODO: use SLIP_matrix_free (&A, option) ;
    // Delete this since *_mat will no longer be used.

    if (A == NULL || *A == NULL) {return;}
    for (int64_t i = m-1; i >= 0; i--)
    {
        SLIP_delete_mpq_array(&((*A)[i]), n);
    }
    SLIP_FREE(*A);
}

