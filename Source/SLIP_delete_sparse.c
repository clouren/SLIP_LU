# include "SLIP_LU_internal.h"

/* Purpose: This function deletes the sparse matrix A */
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
