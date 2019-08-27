# include "SLIP_LU_mex.h"

/* Purpose: This function converts mpq array to double
 * NOTE: This induces roundoff error via the final division
*/
void SLIP_mpq_to_double
(
    double* x_doub,       // double array
    const mpq_t* x_mpq,   // mpq array
    const int32_t n       // size of b
)
{
    SLIP_info status;
    for (int32_t i = 0; i < n; i++)    
    {
        SLIP_MEX_OK(slip_mpq_get_d(&(x_doub[i]), x_mpq[i]));
    }
}
