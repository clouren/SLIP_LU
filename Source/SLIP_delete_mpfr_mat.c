# include "SLIP_LU_internal.h"

/* Purpose: This function deletes a dense mpfr matrix. 
 * 
 * Input is a mpfr*** mat which is destroyed on completion
 */
void SLIP_delete_mpfr_mat
(
    mpfr_t ***A,   // Dense mpfr matrix
    int32_t m,     // number of rows of A
    int32_t n      // number of columns of A
)
{
    if (A == NULL || (*A) == NULL) {return;}
    for (int32_t i = 0; i < m; i++)
    {
        SLIP_delete_mpfr_array(&((*A)[i]), n);
    }
    SLIP_FREE((*A));
}
