#ifndef SLIP_mex 
#define SLIP_mex

#include "SLIP_LU_internal.h"
#include "matrix.h"


#define SLIP_MEX_OK(method)         \
{                                   \
    status = method;                \
    SLIP_mex_error(status);         \
}


/* Purpose: A GMP reallocation function 
 * This allows GMP to use MATLAB's default realloc function 
 */
void* mxGMPRealloc 
(
    void* x,    // void* to be reallocated 
    size_t a,   // Previous size
    size_t b    // New size
);

/* Purpose: A GMP free function. This allows GMP to use
 * MATLAB's mxFree instead of free */

// A GMP realloc function
void mxGMPFree 
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

void SLIP_check_input
(
    const mxArray * input [],    // The matlab input array
    int32_t nargin
);

void SLIP_get_matlab_options
(
    SLIP_options* option,  // Control parameters
    const mxArray* input   // The input options from MATLAB interface
);

void SLIP_int32_to_mwIndex(mwIndex* y, int32_t* x, int32_t n) ;

void SLIP_mwIndex_to_int32(int32_t* y, mwIndex* x, mwSize n) ;

void SLIP_mex_check_for_inf
(
    double* x, // The array of numeric values 
    mwSize n   // size of array
);

/* Purpose: This function reads in the A matrix and right hand side vectors. */
void SLIP_mex_get_A_and_b
(
    SLIP_sparse *A,          // Internal SLIP Mat stored in ccf 
    SLIP_dense *b,           // mpz matrix used internally
    const mxArray* input[],  // The input A matrix and options 
    int32_t nargin           // Number of input to the mexFunction
);


mxArray* SLIP_mex_output_soln(double** x, int32_t m, int32_t n) ;

mxArray* SLIP_mex_output_col_permut(int32_t* x, int32_t n );

mxArray* SLIP_mex_output_p(int32_t* pinv, int32_t n) ;

mxArray* SLIP_mex_output_L
(
    SLIP_sparse *L,    // the sparse matrix to be output
    mpz_t *rhos        // sequence of pivots
);

mxArray* SLIP_mex_output_U
(
    SLIP_sparse *U,    // the sparse matrix to be output
    mpz_t *rhos,       // sequence of pivots
    mpq_t scale        // Scale factor of A matrix

);

void SLIP_mex_error(SLIP_info status) ;

/* This function is modified from CSparse/Source/cs_dropzeros*/
mwIndex SLIP_dropzeros (mxArray *A);

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
mwIndex SLIP_fkeep (mxArray *A, bool (*fkeep) (int32_t, int32_t, double));

/* A = A' */
/* This function is modified from CSparse/Source/cs_transpose*/
SLIP_info SLIP_transpose (mxArray *A);

#endif
