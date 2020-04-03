//------------------------------------------------------------------------------
// SLIP_LU/slip_matrix_div: divide a matrix by a scalar
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function takes as input a dense SLIP_matrix, x, which is MPZ,
 * and divides it by the determinant of the matrix. This division is then stored
 * in a dense MPQ matrix. This is used internally to divide the solution vector 
 * by the determinant of the matrix.
 *
 * On output, the contents of the matrix x2 are modified.
 */

#define SLIP_FREE_ALL             \
    SLIP_MPQ_CLEAR(det2);

#include "slip_LU_internal.h"
SLIP_info slip_matrix_div // divides the x matrix by the determinant
(
    SLIP_matrix* x2,     // solution of x/det
    SLIP_matrix* x,      // input vector
    const mpz_t det      // given determinant of matrix
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    SLIP_REQUIRE (x2, SLIP_DENSE, SLIP_MPQ) ;
    SLIP_REQUIRE (x,  SLIP_DENSE, SLIP_MPZ) ;
    
    //--------------------------------------------------------------------------
    // Set det2 = det
    //--------------------------------------------------------------------------

    mpq_t det2; SLIP_MPQ_SET_NULL(det2);
    SLIP_CHECK(SLIP_mpq_init(det2));
    SLIP_CHECK(SLIP_mpq_set_num(det2, det));

    //--------------------------------------------------------------------------
    // iterate each entry of x, copy to x2 and divide it by det
    //--------------------------------------------------------------------------
    
    int64_t nz = x->n * x->m;
    for (int64_t i = 0; i < nz; i++)
    {
        // Set x2[i] = x[i]
        SLIP_CHECK( SLIP_mpq_set_num( x2->x.mpq[i], x->x.mpz[i]));
        // x2[i] = x2[i] / det2
        SLIP_CHECK(SLIP_mpq_div( x2->x.mpq[i], x2->x.mpq[i], det2));
    }
    
    //--------------------------------------------------------------------------
    // Free memory associated with det2
    //--------------------------------------------------------------------------

    SLIP_FREE_ALL ;
    return SLIP_OK;
}

