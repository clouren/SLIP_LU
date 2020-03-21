//------------------------------------------------------------------------------
// SLIP_LU/SLIP_dense_mpz_mat: delete a dense mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function deletes a dense mpz matrix.
 *
 * Input is a mpz_t*** matrix which is destoyed on completion.
 */

#include "SLIP_LU_internal.h"

void SLIP_delete_mpz_mat
(
    mpz_t ***A,     // The dense mpz matrix
    int32_t m,      // number of rows of A
    int32_t n       // number of columns of A
)
{

    // TODO: use SLIP_matrix_free (&A, option) ;
    // Delete this since *_mat will no longer be used.

    if (A == NULL || (*A) == NULL) {return ;}

    for (int32_t i = 0; i < m; i++)
    {
        SLIP_delete_mpz_array(&((*A)[i]), n);
    }
    SLIP_FREE(*A);
}

