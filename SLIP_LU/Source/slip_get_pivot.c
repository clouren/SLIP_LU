//------------------------------------------------------------------------------
// SLIP_LU/slip_get_pivot: find a pivot entry in a column
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* This function performs the pivoting for the SLIP LU factorization.
 * The optional Order is:
 *
 *  SLIP_SMALLEST = 0,      Smallest pivot
 *  SLIP_DIAGONAL = 1,      Diagonal pivoting
 *  SLIP_FIRST_NONZERO = 2, First nonzero per column chosen as pivot
 *  SLIP_TOL_SMALLEST = 3,  Diagonal pivoting with tolerance for pivot. Default
 *  SLIP_TOL_LARGEST = 4,   Diagonal pivoting with tolerance for largest pivot
 *  SLIP_LARGEST = 5        Largest pivot
 *
 * Options 2 and 5 are not recommended; they are for comparison only.
 *
 * On output, the pivs, rhos, pinv, and row_perm arrays are all modified.
 */

#define SLIP_FREE_ALL           \
    SLIP_MPQ_CLEAR (tol) ;      \
    SLIP_MPQ_CLEAR (ratio) ; 

#include "SLIP_LU_internal.h"

SLIP_info slip_get_pivot
(
    int64_t *pivot, // found index of pivot entry
    mpz_t* x,       // kth column of L and U
    int64_t* pivs,  // vector indicating which rows have been pivotal
    int64_t n,      // dimension of the problem
    int64_t top,    // nonzero pattern is located in xi[top..n-1]
    int64_t* xi,    // nonzero pattern of x
    SLIP_pivot order,  // pivoting method to use (see above description)
    int64_t col,    // current column of A (real kth column i.e., q[k])
    int64_t k,      // iteration of the algorithm
    mpz_t* rhos,    // vector of pivots
    int64_t* pinv,  // row permutation
    int64_t* row_perm,// opposite of pinv. if pinv[i] = j then row_perm[j] = i
    double tolerance// tolerance used if some tolerance based pivoting is used
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    ASSERT (rhos != NULL && rhos->kind == SLIP_DENSE && rhos->type == SLIP_MPZ);

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    int sgn, r;
    mpq_t tol, ratio;
    SLIP_MPQ_SET_NULL(tol);
    SLIP_MPQ_SET_NULL(ratio);

    //--------------------------------------------------------------------------
    // Smallest pivot
    //--------------------------------------------------------------------------

    if (order == SLIP_SMALLEST)
    {
        SLIP_CHECK (slip_get_smallest_pivot(pivot, x, pivs, n, top, xi));
    }

    //--------------------------------------------------------------------------
    // Diagonal
    //--------------------------------------------------------------------------
    else if (order == SLIP_DIAGONAL)
    {
        // Check if x[col] is eligible. take smallest pivot    if not
        SLIP_CHECK (SLIP_mpz_sgn(&sgn, x[col]));
        if (sgn != 0 && pivs[col] < 0)
        {
            *pivot = col;
        }
        else
        {
            SLIP_CHECK (slip_get_smallest_pivot(pivot, x, pivs, n, top, xi));
        }
    }

    //--------------------------------------------------------------------------
    // First nonzero
    //--------------------------------------------------------------------------
    else if (order == SLIP_FIRST_NONZERO)
    {
        SLIP_CHECK (slip_get_nonzero_pivot(pivot, x, pivs, n, top, xi));
    }

    //--------------------------------------------------------------------------
    // Tolerance with largest pivot
    //--------------------------------------------------------------------------
    else if (order == SLIP_TOL_LARGEST)
    {
        SLIP_CHECK (slip_get_largest_pivot(pivot, x, pivs, n, top, xi));

        //----------------------------------------------------------------------
        // Check x[col] vs largest potential pivot
        //----------------------------------------------------------------------
        SLIP_CHECK (SLIP_mpz_sgn(&sgn, x[col]));
        if (sgn != 0 && pivs[col] < 0)
        {
            SLIP_CHECK(SLIP_mpq_init(tol));
            SLIP_CHECK(SLIP_mpq_init(ratio));
            // tol = user specified tolerance
            SLIP_CHECK(SLIP_mpq_set_d(tol, tolerance));
            // ratio = diagonal/largest
            SLIP_CHECK(SLIP_mpq_set_num(ratio, x[col]));
            SLIP_CHECK(SLIP_mpq_set_den(ratio, x[*pivot]));
            // ratio = |ratio|
            SLIP_CHECK(SLIP_mpq_abs(ratio, ratio));

            // Is ratio >= tol?
            SLIP_CHECK(SLIP_mpq_cmp(&r, ratio, tol));
            if (r >= 0)
            {
                *pivot = col;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Use the largest potential pivot
    //--------------------------------------------------------------------------
    else if (order == SLIP_LARGEST)
    {
        SLIP_CHECK (slip_get_largest_pivot(pivot, x, pivs, n, top, xi));
    }

    //--------------------------------------------------------------------------
    // Tolerance with smallest pivot (default option)
    //--------------------------------------------------------------------------
    else // if (order == SLIP_TOL_SMALLEST)
    {
        SLIP_CHECK (slip_get_smallest_pivot(pivot, x, pivs, n, top, xi)) ;

        //----------------------------------------------------------------------
        // Checking x[col] vs smallest pivot
        //----------------------------------------------------------------------
        SLIP_CHECK (SLIP_mpz_sgn(&sgn, x[col]));
        if (sgn != 0 && pivs[col] < 0)
        {

            // Initialize tolerance and ratio
            SLIP_CHECK(SLIP_mpq_init(tol));
            SLIP_CHECK(SLIP_mpq_init(ratio));

            // ratio = |smallest/diagonal|
            SLIP_CHECK(SLIP_mpz_abs(SLIP_MPQ_NUM(ratio), x[*pivot]));
            SLIP_CHECK(SLIP_mpz_abs(SLIP_MPQ_DEN(ratio), x[col]));

            // Set user specified tolerance
            SLIP_CHECK(SLIP_mpq_set_d(tol, tolerance));

            // Is ratio >= tol?
            SLIP_CHECK(SLIP_mpq_cmp(&r, ratio, tol));
            if (r >= 0)
            {
                *pivot = col;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Reflect changes in row location & row_perm
    //--------------------------------------------------------------------------
    // Must move pivot into position k
    int64_t intermed = pinv[*pivot];
    int64_t intermed2 = row_perm[k];

    //--------------------------------------------------------------------------
    // Set row_perm[k] = pivot and row_perm[pinv[pivot]] = row_perm[k]
    // Also, set pinv[pivot] = k and pinv[row_perm[k]] = pinv[pivot]
    //--------------------------------------------------------------------------
    row_perm[k] = *pivot;
    row_perm[intermed] = intermed2;
    pinv[*pivot] = k;
    pinv[intermed2] = intermed;
    // Row pivot is now pivotal
    pivs[*pivot] = 1;
    // The kth pivot is x[pivot]
    SLIP_CHECK (SLIP_mpz_set(rhos[k], x[*pivot]));

    // Free memory
    SLIP_FREE_ALL;
    return SLIP_OK;
}

