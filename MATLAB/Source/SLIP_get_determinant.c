# include "SLIP_LU_mex.h"

/* This function gets the determinant of the input matrix
 * It is utilized within the SLIP LU matlab interface
 * 
 * on output the mpq_t deter contains the determinant
 */

void SLIP_get_determinant
(
    mpq_t deter,            // output is mpq_t determinant
    mpz_t det,              // determinant of scaled A
    SLIP_sparse *A          // sparse matrix A
)
{
    if (!A)
    {
    	SLIP_mex_error (SLIP_INCORRECT_INPUT);
    }
    int32_t r;
    SLIP_info status;

    // Set deter = deter of scaled A
    SLIP_MEX_OK(slip_mpq_set_z(deter, det));

    // Is LU scale 1?
    SLIP_MEX_OK(slip_mpq_cmp_ui(&r, A->scale, 1, 1)); 
    if (r != 0)
    {
        // Scale deter
        for (int32_t k = 0; k < A->n; k++)
        {
            SLIP_MEX_OK(slip_mpq_div(deter, deter, A->scale));
        }
    }
}
