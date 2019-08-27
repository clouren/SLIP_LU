# include "SLIP_LU_internal.h"

/* Purpose: This function frees the memory of the SLIP_col struct
 *
 * Input is the SLIP_col structure, it is destroyed on function termination.
 */
void SLIP_delete_LU_analysis
(
    SLIP_LU_analysis **S // Structure to be deleted 
)
{
    if ((S == NULL) || (*S == NULL)) {return;}
    SLIP_FREE((*S)->q);
    SLIP_FREE(*S);
}
