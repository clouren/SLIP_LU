#include "SLIP_LU_mex.h"


/* Purpose: This function outputs a double array as a mxArray. */
mxArray* SLIP_mex_output_determinant
(
    mpz_t scaled_det,       // determinant of scaled A
    SLIP_sparse *A          // sparse matrix A
)
{
    SLIP_info status;
    // Create a n * 1 array
    mxArray* Xmatlab = mxCreateDoubleMatrix (1, 1, mxREAL);

    // Get the numeric values
    double* d = mxGetPr(Xmatlab);

    mpq_t deter;
    SLIP_MPQ_SET_NULL(deter);
    SLIP_MEX_OK (slip_mpq_init(deter));
    SLIP_get_determinant(deter, scaled_det, A);
    SLIP_MEX_OK (slip_mpq_get_d(&(d[0]), deter));

    d[0] = fabs(d[0]);
    SLIP_MPQ_CLEAR(deter);

    return Xmatlab;
}
