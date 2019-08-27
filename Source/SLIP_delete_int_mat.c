# include "SLIP_LU_internal.h"

/* Purpose: This function deletes a dense int matrix.
 * 
 * Input is a int*** mat and its dimensions.
 * Input mat is destroyed on function completion
 */
void SLIP_delete_int_mat
(
    int32_t ***A,  // dense matrix
    int32_t m,     // number of rows
    int32_t n      // number of columns
)
{
    if (A == NULL || (*A) == NULL) {return;}
    for (int32_t i = 0; i < m; i++)
    {
        SLIP_FREE( (*A)[i]);
    }
    SLIP_FREE(*A);
}
