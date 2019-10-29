#include "SLIP_LU_mex.h"

/* Purpose: This function outputs the sparse matrix U*/
mxArray* SLIP_mex_output_U
(
    SLIP_sparse *U,    // the sparse matrix to be output
    mpz_t *rhos,       // sequence of pivots
    mpq_t scale        // Scale factor of A matrix
)
{
    if (!U)  { SLIP_mex_error (SLIP_INCORRECT_INPUT); }
    SLIP_info status;
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_MEX_OK(slip_mpq_init(temp));

    // Create a m*n sparse matrix
    mxArray *Amatlab = mxCreateSparse((mwSize) U->m, (mwSize) U->n,
                                      (mwSize) U->nz, mxREAL);
    mwIndex *pA, *iA;
    double *xA;
    pA = mxGetJc(Amatlab);
    iA = mxGetIr(Amatlab);
    xA = mxGetDoubles(Amatlab);

    // Populate U->p and U->i
    SLIP_int32_to_mwIndex(pA, U->p, (U->n)+1);
    SLIP_int32_to_mwIndex(iA, U->i, U->nz);

    //--------------------------------------------------------------------------
    // As per Escobedo et. al., row j of U is divided by rho[j-1]
    //--------------------------------------------------------------------------
    mpq_t* xU = SLIP_create_mpq_array(U->nz);
    if (!xU)  { SLIP_MEX_OK (SLIP_OUT_OF_MEMORY); }

    for (int32_t j = 0; j < U->nz; j++)
    {
        SLIP_MEX_OK(slip_mpq_set_z(xU[j], U->x[j]));
        if (U->i[j] > 0)
        {
            SLIP_MEX_OK(slip_mpq_set_z(temp, rhos[U->i[j]-1]));             
            SLIP_MEX_OK(slip_mpq_div(xU[j], xU[j], temp));
        }
        // Populate xU
        SLIP_MEX_OK(slip_mpq_div(xU[j], xU[j], scale));
    }

    SLIP_mpq_to_double(xA, xU, U->nz);
    SLIP_delete_mpq_array(&xU, U->nz);
    SLIP_MPQ_CLEAR(temp);

    // drop zeros from U and sort it
    SLIP_dropzeros(Amatlab);
    SLIP_transpose(Amatlab);
    SLIP_transpose(Amatlab);
    return Amatlab;
}
