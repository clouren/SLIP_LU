#include "SLIP_LU_mex.h"

/* Purpose: This function outputs the sparse matrix L*/
mxArray* SLIP_mex_output_L
(
    SLIP_sparse *L,    // the sparse matrix to be output
    mpz_t *rhos        // sequence of pivots
)
{
    if (!L)  { SLIP_mex_error (SLIP_INCORRECT_INPUT); }
    SLIP_info status;
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_MEX_OK(slip_mpq_init(temp));

    // Create a m*n sparse matrix
    mxArray *Amatlab = mxCreateSparse((mwSize) L->m, (mwSize) L->n,
                                      (mwSize) L->nz, mxREAL);
    mwIndex *pA, *iA;
    double *xA;
    pA = mxGetJc(Amatlab);
    iA = mxGetIr(Amatlab);
    xA = mxGetDoubles(Amatlab);

    // Populate L->p and L->i
    SLIP_int32_to_mwIndex(pA, L->p, (L->n)+1);
    SLIP_int32_to_mwIndex(iA, L->i, L->nz);

    //--------------------------------------------------------------------------
    // As per Escobedo et. al. Column j of L is divided by rho[j]
    //--------------------------------------------------------------------------
    mpq_t* xL = SLIP_create_mpq_array(L->nz);
    if (!xL)  { SLIP_MEX_OK (SLIP_OUT_OF_MEMORY); }

    for (int32_t k = 0; k < L->n; k++)
    {
        // Populate pL
        SLIP_MEX_OK(slip_mpq_set_z(temp, rhos[k]));
        for (int32_t j = L->p[k]; j < L->p[k+1]; j++)
        {
            SLIP_MEX_OK(slip_mpq_set_z(xL[j], L->x[j]));
            // Populate xL
            SLIP_MEX_OK(slip_mpq_div(xL[j], xL[j], temp));
        }
    }

    SLIP_mpq_to_double(xA, xL, L->nz);
    SLIP_delete_mpq_array(&xL, L->nz);
    SLIP_MPQ_CLEAR(temp);

    // drop zeros from L and sort it
    SLIP_dropzeros(Amatlab);
    SLIP_transpose(Amatlab);
    SLIP_transpose(Amatlab);
    return Amatlab;
}
