//------------------------------------------------------------------------------
// SLIP_LU/slip_scale_matrix: scale a matrix if its coming from mpz values
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// When making a copy of matrix, if the output matrix is either mpq, mpfr,
// int64_t, or fp64, the scaling factor associated with this matrix must be 1.
// Scaling factors are only valid for mpz_t matrices. This function serves as
// the final component of SLIP_matrix_copy. It accepts as input the final
// (copied) C matrix and the original A matrix. It then applies the scaling
// factor if necessary.

#define SLIP_FREE_ALL   \
SLIP_MPQ_CLEAR(temp);   \

#include "SLIP_LU_internal.h"

SLIP_info slip_scale_matrix
(
    SLIP_matrix *C,             // Copied (output) matrix
    SLIP_matrix *A,             // Matrix which was copied (input)
    SLIP_options *option        // Command options, if NULL defaults are used
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info = SLIP_OK ;
    
    if (C->type == SLIP_MPZ || // Already has been accounted for
        A->type != SLIP_MPZ || // A is not mpz means A has no scaling factor
        C->x_shallow)          // Not modifying shallow arrays Matrix copy is
                               // never shallow So this should never be
                               // triggered TODO but C is never shallow
    {
        return info;
    }
    
    int r;
    int64_t i;
    
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(SLIP_mpq_init(temp));

    // Determine if A's scaling factor is 1
    SLIP_CHECK(SLIP_mpq_cmp_ui(&r, A->scale, 1, 1));
    if (r == 0)     // No scaling factor applied
        return info;
    
    switch (C->type)
    {
        case SLIP_MPQ:
        {
            // C is MPQ, scale is mpq, do a simple divide.
            for (i = 0; i < C->nzmax; i++)
            {
                SLIP_CHECK(SLIP_mpq_div( C->x.mpq[i], C->x.mpq[i], A->scale));
            }
            break;
        }
        case SLIP_MPFR:
        {
            // C is MPFR, scale is mpq. We will cast each mpfr to mpq, perform
            // the division, then cast back to mpfr.
            
            for (i = 0; i < C->nzmax; i++)
            {
                // temp = C->x[i]
                SLIP_CHECK(SLIP_mpfr_get_q(temp, C->x.mpfr[i], option->round));
                // temp = temp/scale
                SLIP_CHECK(SLIP_mpq_div(temp, temp, A->scale));
                // C->x[i] = temp
                SLIP_CHECK(SLIP_mpfr_set_q(C->x.mpfr[i], temp, option->round));
            }
            break;
        }
        case SLIP_INT64:
        {
            // C is INT64, scale is mpq. We will cast each int64 to mpq, perform
            // the division, then cast back to int64_t.
            // WARNING: This leads to truncation and loss of data.
            for (i = 0; i < C->nzmax; i++)
            {
                // temp = C->x[i]
                SLIP_CHECK (SLIP_mpq_set_si (temp, C->x.int64 [i], 1)) ;
                // temp = temp/scale
                SLIP_CHECK(SLIP_mpq_div(temp, temp, A->scale));
                // C->x[i] = temp
                double t;
                SLIP_CHECK (SLIP_mpq_get_d (&t, temp)) ;
                C->x.int64[i] = slip_cast_double_to_int64 (t) ;
            }
            break;
        }
        case SLIP_FP64:
        {
            // C is double, scale is mpq. We will cast each double to mpq,
            // perform the division, then cast back to double.
            
            for (i = 0; i < C->nzmax; i++)
            {
                // temp = C->x[i]
                SLIP_CHECK(SLIP_mpq_set_d(temp, C->x.fp64[i]));
                // temp = temp/scale
                SLIP_CHECK(SLIP_mpq_div(temp, temp, A->scale));
                // C->x[i] = temp
                SLIP_CHECK(SLIP_mpq_get_d( &C->x.fp64[i], temp));
            }
            break;
        }
        default: SLIP_FREE_ALL; return SLIP_INCORRECT_INPUT;
    }    

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    SLIP_mpq_set_ui(C->scale, 1, 1);
            
    SLIP_CHECK (info) ;
    SLIP_FREE_ALL;
    return (SLIP_OK) ;
}

