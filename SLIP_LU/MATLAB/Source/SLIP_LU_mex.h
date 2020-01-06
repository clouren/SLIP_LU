//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/SLIP_LU_mex.h: include file for MATLAB functions
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#ifndef SLIP_mex 
#define SLIP_mex

#include "SLIP_LU_internal.h"
#include "matrix.h"


#define SLIP_MEX_OK(method)         \
{                                   \
    status = method;                \
    slip_mex_error(status);         \
}


/* Purpose: A GMP reallocation function 
 * This allows GMP to use MATLAB's default realloc function 
 */
void* SLIP_gmp_mex_realloc 
(
    void* x,    // void* to be reallocated 
    size_t a,   // Previous size
    size_t b    // New size
);

/* Purpose: A GMP free function. This allows GMP to use
 * MATLAB's mxFree instead of free 
 */
void SLIP_gmp_mex_free 
(
    void* x,    // void* to be freed
    size_t a    // Size
);

/* Purpose: This function converts mpq array to double
 * NOTE: This induces roundoff error via the final division
*/
void SLIP_mpq_to_double
(
    double* x_doub,       // double array
    const mpq_t* x_mpq,   // mpq array
    const int32_t n       // size of b
);

void slip_check_input
(
    const mxArray * input [],    // The matlab input array
    int32_t nargin
);

void slip_get_matlab_options
(
    SLIP_options* option,  // Control parameters
    const mxArray* input   // The input options from MATLAB interface
);

/* Purpose: Convert int32_t* array to MATLAB mwIndex* array
 */
void slip_int32_to_mwIndex
(
    mwIndex* y, 
    int32_t* x, 
    int32_t n
) ;

/* Purpose: convert mwIndex* array to int32_t* array
 */
void slip_mwIndex_to_int32
(
    int32_t* y, 
    mwIndex* x, 
    mwSize n
) ;

/* Purpose: Check the input x array for numbers too large for 
 * double precision.
 */
void slip_mex_check_for_inf
(
    double* x, // The array of numeric values 
    mwSize n   // size of array
);

/* Purpose: This function reads in the A matrix and right hand side vectors. */
void slip_mex_get_A_and_b
(
    SLIP_sparse *A,          // Internal SLIP Mat stored in ccf 
    SLIP_dense *b,           // mpz matrix used internally
    const mxArray* input[],  // The input A matrix and options 
    int32_t nargin           // Number of input to the mexFunction
);


/* Purpose: Output the solution to the linear system Ax=b to matlab
 */
mxArray* slip_mex_output_soln
(
    double** x, 
    int32_t m, 
    int32_t n
) ;

/* Purpose: Output the column permutation the Q matrix 
 */
mxArray* slip_mex_output_col_permut
(
    int32_t* x, 
    int32_t n 
);

/* Purpose: Output the row permutation the P matrix 
 */
mxArray* slip_mex_output_p
(
    int32_t* pinv, 
    int32_t n
) ;

/* Purpose: Output the L matrix 
 */
mxArray* slip_mex_output_L
(
    SLIP_sparse *L,    // the sparse matrix to be output
    mpz_t *rhos        // sequence of pivots
);

/* Purpose: Output the U matrix 
 */
mxArray* slip_mex_output_U
(
    SLIP_sparse *U,    // the sparse matrix to be output
    mpz_t *rhos,       // sequence of pivots
    mpq_t scale        // Scale factor of A matrix

);

/* Purpose: Report errors if they arise
 */
void slip_mex_error
(
    SLIP_info status
) ;

/* Purpose: Drop entries which are zero from a sparse matrix
 */
mwIndex slip_dropzeros 
(
    mxArray *A
);

/* Purpose: Used to drop zeros 
 */
mwIndex slip_fkeep 
(
    mxArray *A, 
    bool (*fkeep) (int32_t, int32_t, double)
);

/* Purpose: transpose the matrix A
 */
SLIP_info slip_transpose 
(
    mxArray *A
);

#endif
