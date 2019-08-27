# include "SLIP_LU_internal.h"

/* Purpose: This function deletes a dense mpq matrix
 *
 * Input is a mpq_t*** matrix which is destroyed upon function completion
 */
void SLIP_delete_mpq_mat
(
    mpq_t***A,     // dense mpq matrix
    int32_t m,     // number of rows of A
    int32_t n      // number of columns of A
)
{
    //--------------------------------------------------------------------------
    // Iterate accross all entries and clear the individual memory
    //--------------------------------------------------------------------------
    if (A == NULL || *A == NULL) {return;}
    for (int32_t i = m-1; i >= 0; i--)
    {
        SLIP_delete_mpq_array(&((*A)[i]), n);
    }
    SLIP_FREE(*A);
}
