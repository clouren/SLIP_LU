//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_sparse: delete a sparse mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function deletes the sparse matrix A */

#include "SLIP_LU_internal.h"

void SLIP_delete_sparse
(
    SLIP_sparse ** A // matrix to be deleted
)
{
    if (A == NULL || (*A) == NULL) {return ;}
    SLIP_delete_mpz_array(&((*A)->x), (*A)->nzmax);
    SLIP_FREE ((*A)->i);
    SLIP_FREE ((*A)->p);
    SLIP_MPQ_CLEAR( (*A)->scale);
    SLIP_FREE ((*A));
}

