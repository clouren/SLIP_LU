//------------------------------------------------------------------------------
// SLIP_LU/SLIP_delete_dense: destroy a dense mpz matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: Delete the SLIP_dense data structure */

void SLIP_delete_dense
( 
    SLIP_dense **A   // Struct to be destroyed
)
{
    if (A == NULL || (*A) == NULL){ return;}
 
    // Delete A->x
    SLIP_delete_mpz_mat(&((*A)->x), (*A)->m, (*A)->n);
    // Delete A->scale
    SLIP_MPQ_CLEAR((*A)->scale);
    SLIP_FREE(*A);
}
