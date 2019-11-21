#ifndef SLIP_Include
#define SLIP_Include

// This software package exactly solves a sparse system of linear equations
// using the SLIP LU factorization. This code accompanies the paper (submitted
// to ACM Transactions on Mathematical Software):

//    "Algorithm 1xxx: SLIP LU: A Sparse Left-Looking Integer-Preserving LU
//    Factorization for Exactly Solving Sparse Linear Systems",
//    C. Lourenco, J. Chen, E. Moreno-Centeno, T. Davis, under submission,
//    ACM Trans. Mathematical Software.

//    The theory associated with this software can be found in the paper
//    (published in SIAM journal on matrix analysis and applications):

//    "Exact Solution of Sparse Linear Systems via Left-Looking
//     Roundoff-Error-Free LU Factorization in Time Proportional to
//     Arithmetic Work", C. Lourenco, A. R. Escobedo, E. Moreno-Centeno,
//     T. Davis, SIAM J. Matrix Analysis and Applications.  pp 609-638,
//     vol 40, no 2, 2019.

//    If you use this code, you must first download and install the GMP and
//    MPFR libraries. GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/

//    If you use SLIP LU for a publication, we request that you please cite
//    the above two papers.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco, Jinhao Chen, Erick Moreno-Centeno, and Timothy Davis

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SLIP LU is free software; you can redistribute it and/or modify
//     it under the terms of either:
//
//        * the GNU Lesser General Public License as published by the
//          Free Software Foundation; either version 3 of the License,
//          or (at your option) any later version.
//
//     or
//
//        * the GNU General Public License as published by the Free Software
//          Foundation; either version 2 of the License, or (at your option) any
//          later version.
//
//    or both in parallel, as here.
//
//    See license.txt for license info.
//
// This software is copyright by Christopher Lourenco, Jinhao Chen, Erick
// Moreno-Centeno and Timothy A. Davis. All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SLIP LU is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    This software package solves the linear system Ax = b exactly. The input
//    matrix and right hand side vectors are stored as either integers, double
//    precision numbers, multiple precision floating points (through the mpfr
//    library) or as rational numbers (as a collection of numerators and
//    denominators using the GMP mpq_t data structure). Appropriate routines
//    within the code transform the input into an integral matrix in compressed
//    column form.

//    This package computes the factorization PAQ = LDU. Note that we store the
//    "functional" form of the factorization by only storing L and U. The user
//    is given some freedom to select the permutation matrices P and Q. The
//    recommended default settings select Q using the COLAMD column ordering
//    and select P via a partial pivoting scheme in which the diagonal entry
//    in column k is selected if it is the same magnitude as the smallest
//    entry, otherwise the smallest entry is selected as the kth pivot.
//    Alternative strategies allowed to select Q include the AMD column
//    ordering or no column permutation (Q=I).  For pivots, there are a variety
//    of potential schemes including traditional partial pivoting, diagonal
//    pivoting, tolerance pivoting etc. This package does not allow pivoting
//    based on sparsity criterion.

//    The factors L and U are computed via integer preserving operations via
//    integer-preserving Gaussian elimination. The key part of this algorithm
//    is a REF Sparse triangular solve function which exploits sparsity to
//    reduce the number of operations that must be performed.

//    Once L and U are computed, a simplified version of the triangular solve
//    is performed which assumes the vector b is dense. The final solution
//    vector x is gauranteed to be exact. This vector can be output in one of
//    three ways: 1) full precision rational arithmetic (as a sequence of
//    numerators and denominators) using the GMP mpq_t data type, 2) double
//    precision while not exact will produce a solution accurate to machine
//    roundoff unless the size of the associated solution exceeds double
//    precision (i.e., the solution is 10^500 or something), 3) variable
//    precision floating point using the GMP mpfr_t data type. The associated
//    precision is user defined.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Current version of the code
#define SLIP_LU_VERSION "1.0.0"
#define SLIP_LU_VERSION_MAJOR 1
#define SLIP_LU_VERSION_MINOR 0
#define SLIP_LU_VERSION_SUB   0

// Name of associated paper
#define SLIP_PAPER "Algorithm 1xxx: SLIP LU: Sparse Left-looking Integer-Preserving LU Factorization"

// Authors of code
#define SLIP_AUTHOR "Christopher Lourenco, Jinhao Chen, Erick Moreno-Centeno, Timothy Davis"

//------------------------------------------------------------------------------
// Type of MPFR rounding used. 
//------------------------------------------------------------------------------

// The MPFR library utilizes an internal rounding scheme. The options are
//  MPFR_RNDN: round to nearest (roundTiesToEven in IEEE 754-2008),
//  MPFR_RNDZ: round toward zero (roundTowardZero in IEEE 754-2008),
//  MPFR_RNDU: round toward plus infinity (roundTowardPositive in IEEE 754-2008),
//  MPFR_RNDD: round toward minus infinity (roundTowardNegative in IEEE 754-2008),
//  MPFR_RNDA: round away from zero.
//  MPFR_RNDF: faithful rounding. This is not stable

// SLIP LU utilizes MPFR_RNDN. If the user wishes to change this, they can change 
// the following parameter

#define SLIP_MPFR_ROUND MPFR_RNDN


//------------------------------------------------------------------------------
// Error codes
//------------------------------------------------------------------------------

// ALL SLIP_LU functions return a code that indicates if it was successful
// or not.

typedef enum
{
    SLIP_OK = 0,                    // all is well
    SLIP_OUT_OF_MEMORY = -1,        // out of memory
    SLIP_SINGULAR = -2,             // the input matrix A is singular
    SLIP_INCORRECT_INPUT = -3,      // one or more input arguments are incorrect
    SLIP_INCORRECT = -4             // The solution is incorrect
}
SLIP_info ;

//------------------------------------------------------------------------------
// Pivot scheme codes
//------------------------------------------------------------------------------

// A code in SLIP_options to tell SLIP LU what type of pivoting to use.

typedef enum
{
    SLIP_SMALLEST = 0,              // Smallest pivot
    SLIP_DIAGONAL = 1,              // Diagonal pivoting
    SLIP_FIRST_NONZERO = 2,         // First nonzero per column chosen as pivot
    SLIP_TOL_SMALLEST = 3,          // Diagonal pivoting with tolerance for
                                    // smallest pivot. Default
    SLIP_TOL_LARGEST = 4,           // Diagonal pivoting with tolerance for
                                    // largest pivot
    SLIP_LARGEST = 5                // Largest pivot
}
SLIP_pivot ;

//------------------------------------------------------------------------------
// Column ordering scheme codes
//------------------------------------------------------------------------------

// A code in SLIP_options to tell SLIP LU what column ordering to use.

typedef enum
{
    SLIP_NO_ORDERING = 0,           // None: Not recommended for sparse matrices
    SLIP_COLAMD = 1,                // COLAMD: Default
    SLIP_AMD = 2                    // AMD
}
SLIP_col_order ;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Data Structures--------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/* This struct defines the command line options for the factorization.
 * In addition, it stores the determinant of the matrix and
 * scaling parameter used on A and b. In essence it serves as a global struct to
 * define all options
 */

typedef struct SLIP_options
{
    bool check;         // TRUE if the solution will be checked
    SLIP_pivot pivot;   // Type of pivoting scheme used.
    SLIP_col_order order;// Type of column ordering scheme used
    double tol;         // User specified tolerance for SLIP_TOL_SMALLEST and
                        // SLIP_TOL_LARGEST
    int32_t print_level;// 0: print nothing, 1: just errors,
                        // 2: terse, with basic stats from COLAMD/AMD and SLIP
                        // 3: all, with matrices and results
    uint64_t prec;      // Precision used to output file if MPFR is chosen
} SLIP_options;

/* Purpose: Create and return SLIP_options pointer with default parameters
 * upon successful allocation, which are defined in SLIP_LU_internal.h
 * To free it, simply use SLIP_FREE(option)
 */
SLIP_options* SLIP_create_default_options(void);

//------------------------------------------------------------------------------
// SLIP_sparse: a sparse matrix in compressed sparse column form
//------------------------------------------------------------------------------

// This struct defines a matrix stored in sparse compressed column form. Since
// this code deals with exact data structures, each SLIP_sparse stores
// both the size of internal vectors and the number of elements allocated in
// each one.  This is done to use as little memory as possible.

// The matrix in sparse compressed column form is stored as follows.  p is an
// array of size n+1.  The row indices of entries in column j are held in A->i
// [A->p [j] ... A->p [j+1]-1].  The corresponding values are in the same
// locations in A->x.  A->scale defines how the original input was scaled to
// create the mpz_t integer form.  For an L or U factor, the scale is 1.

typedef struct
{
    int32_t m;    // Number of rows. Current version always uses n == m
    int32_t n;    // Number of columns
    int32_t nzmax;// Allocated size of A->i and A->x
    int32_t nz;   // Number of nonzeros in the matrix
    int32_t *p;   // Column pointers. Array size is n+1
    int32_t *i;   // Row indices. Array size is nzmax, # of entries = nz
    mpz_t *x;     // Values in matrix with size of nzmax, # of entries = nz
    mpq_t scale ; // scale factor for the matrix
} SLIP_sparse ;

/* Purpose: This function return an created empty sparse matrix as SLIP_sparse
 * pointer upon successful malloc
 */
SLIP_sparse *SLIP_create_sparse( void );

/* Purpose: This function deletes the sparse matrix A */
void SLIP_delete_sparse
(
    SLIP_sparse **A // matrix to be deleted
);

//------------------------------------------------------------------------------
// SLIP_dense: a dense 2D matrix of mpz_t entries
//------------------------------------------------------------------------------

// If B is pointer to a m-by-n SLIP_dense matrix, then B->x [i][j] is the
// entry B(i,j), of type mpz_t. The scaling factor is how the initial input
// was modified in order to make B integral.

typedef struct
{
    int32_t m;    // Number of rows
    int32_t n;    // Number of columns
    mpz_t **x;    // Values in matrix with size of m*n
    mpq_t scale ; // scale factor for the matrix

} SLIP_dense ;

/* Purpose: This function creates an empty SLIP_dense matrix */
SLIP_dense *SLIP_create_dense( void );

/* Purpose: This function deletes the dense matrix A */
void SLIP_delete_dense
(
    SLIP_dense **A
);


//------------------------------------------------------------------------------
// SLIP_LU_analysis: symbolic pre-analysis
//------------------------------------------------------------------------------

/* This struct stores the column permutation for LU and the guess on nnz for
 * L and U */

typedef struct
{
    int32_t *q;     // Column permutation for LU
    int32_t lnz;    // Approximate number of nonzeros in L.
                    // i.e., initial size of L
    int32_t unz;    // Approximate number of nonzeros in U.
                    // i.e., initial size of U
} SLIP_LU_analysis;

/* Purpose: This function returns a pointer to a created SLIP_LU_analysis type
 * with the length of S->q set as n (which needs to 1 + number of rows of input
 * matrix) upon successful malloc, otherwise, return NULL
 */
SLIP_LU_analysis *SLIP_create_LU_analysis
(
    int32_t n      // length of S->q
);

/* Purpose: This function frees the memory of the SLIP_LU_analysis struct
 *
 * Input is the SLIP_LU_analysis structure, it is destroyed on function
 * termination.
 */
void SLIP_delete_LU_analysis
(
    SLIP_LU_analysis **S // Structure to be deleted
);

//------------------------------------------------------------------------------
// Memory management
//------------------------------------------------------------------------------

/*
 * Purpose: calloc space of size n*size
 * on failure, NULL is returned
 */

void * SLIP_calloc
(
    size_t n,          // Size of array
    size_t size        // Size to alloc
);

/* Purpose: Define malloc and free for SLIP LU
 *
 * Output arguments are not modified, returned is either a pointer to
 * size space or a NULL pointer in the case of failure
 */
void * SLIP_malloc
(
    size_t size        // Size to alloc
);

/* If p is non-NULL on input, it points to a previously allocated object of
 * size old_size * size_of_item.  The object is reallocated to be of size
 * new_size * size_of_item.  If p is NULL on input, then a new object of that
 * size is allocated.  On success, a pointer to the new object is returned, and
 * ok is returned as true.  If the allocation fails, ok is set to false and a
 * pointer to the old (unmodified) object is returned.
*/
void* SLIP_realloc
(
    void *p,            // Pointer to be realloced
    size_t old_size,    // Old size of this pointer
    size_t new_size     // New size of this pointer
);

/* Purpose: Free the memory associated with the pointer x
 *
 * If we have defined MATLAB, then we use MATLAB's mxFree,
 * otherwise, we use default free.
 *
 * This function safely does nothing if x is NULL on input.
 */
void SLIP_free
(
    void *p         // Pointer to be free'd
);

// Free a pointer and set it to NULL.
#define SLIP_FREE(p)                        \
{                                           \
    SLIP_free (p) ;                         \
    (p) = NULL ;                            \
}

//------------------------------------------------------------------------------
//---------------------------Build the users input matrix from ccf--------------
//------------------------------------------------------------------------------

/* SLIP_build_sparse_ccf_[type] will allow the user to take a matrix of their
 * defined type (either double, mpfr_t, mpz_t, or mpq_t) and convert it from
 * their version of compressed column form to our data structure. The integrity
 * of the user defined arrays are maintained (therefore, one would need to
 * delete these arrays)
 *
 * On output, the SLIP_sparse* A structure contains the input matrix
 *
 */

/* Purpose: Build a SLIP_sparse from mpz_t stored ccf matrix */

SLIP_info SLIP_build_sparse_ccf_mpz
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    mpz_t *x,             // Set of values in full precision int.
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
);

/* Purpose: Build a SLIP_sparse from double stored ccf matrix */

SLIP_info SLIP_build_sparse_ccf_double
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    double *x,            // Set of values as doubles
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
);

/* Purpose: Build a SLIP_sparse from int stored ccf matrix */

SLIP_info SLIP_build_sparse_ccf_int
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    int32_t *x,           // Set of values as doubles
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
);

/* Purpose: Build a SLIP_sparse from mpq_t stored ccf matrix */

SLIP_info SLIP_build_sparse_ccf_mpq
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    mpq_t *x,             // Set of values as mpq_t rational numbers
    int32_t n,            // dimension of the matrix
    int32_t nz            // number of nonzeros in A (size of x and I vectors)
);

/* Purpose: Build a SLIP_sparse from int stored ccf matrix */

SLIP_info SLIP_build_sparse_ccf_mpfr
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *p,           // The set of column pointers
    int32_t *I,           // set of row indices
    mpfr_t *x,            // Set of values as doubles
    int32_t n,            // dimension of the matrix
    int32_t nz,           // number of nonzeros in A (size of x and I vectors)
    SLIP_options *option  // command options containing the prec for mpfr
);

//------------------------------------------------------------------------------
//---------------------------Build the users input matrix from triplet----------
//------------------------------------------------------------------------------
/* SLIP_build_sparse_trip_[type] will allow the user to take a matrix of their
 * defined type (either double, mpfr_t, mpz_t, or mpq_t) and convert it from
 * their triplet form to our data structure. The integrity of the user defined
 * arrays are maintained (therefore, one would need to delete these arrays)
 *
 * On output, the SLIP_sparse* A contains the user's matrix
 *
 */

/* Purpose: Build a sparse matrix from triplet form where input is mpz_t */
SLIP_info SLIP_build_sparse_trip_mpz
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    mpz_t *x,           // Set of values in full precision int
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x, I, and J vectors)
);

/* Purpose: Build a sparse matrix from triplet form where input is double */
SLIP_info SLIP_build_sparse_trip_double
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    double *x,          // Set of values in double
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x, I, and J vectors)
);

/* Purpose: Build a sparse matrix from triplet form where input is int */
SLIP_info SLIP_build_sparse_trip_int
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    int32_t *x,             // Set of values in int
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x, I, and J vectors)
);

/* Purpose: Build a sparse matrix from triplet form where input is mpq_t */
SLIP_info SLIP_build_sparse_trip_mpq
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    mpq_t *x,           // Set of values as rational numbers
    int32_t n,          // dimension of the matrix
    int32_t nz          // number of nonzeros in A (size of x, I, and J vectors)
);

/* Purpose: Build a sparse matrix from triplet form where input is mpfr_t */
SLIP_info SLIP_build_sparse_trip_mpfr
(
    SLIP_sparse *A_output,// It should be initialized but unused yet
    int32_t *I,         // set of row indices
    int32_t *J,         // set of column indices
    mpfr_t *x,          // Set of values as mpfr_t
    int32_t n,          // dimension of the matrix
    int32_t nz,         // number of nonzeros in A (size of x, I, and J vectors)
    SLIP_options *option// command options containing the prec for mpfr
);

//------------------------------------------------------------------------------
//---------------------------Build a dense matrix-------------------------------
//------------------------------------------------------------------------------

/* SLIP_build_dense_[type] will allow the user to take a rhs matrix of their
 * defined type (either double, mpfr_t, mpz_t, or mpq_t) and convert it from
 * their form to mpz. The integrity of the user defined arrays are maintained
 * (therefore, one would need to delete these arrays).
 *
 * On output, the SLIP_dense *A contains the user's matrix
 *
 */

/* Purpose: Build a dense matrix from mpz input */
SLIP_info SLIP_build_dense_mpz
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    mpz_t **b,            // Set of values in full precision int.
    int32_t m,            // number of rows
    int32_t n             // number of columns
);

/* Purpose: Build a dense matrix from double input */
SLIP_info SLIP_build_dense_double
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    double **b,           // Set of values as doubles
    int32_t m,            // number of rows
    int32_t n             // number of columns
);

/* Purpose: Build a dense matrix from int input */
SLIP_info SLIP_build_dense_int
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    int32_t **b,          // Set of values as ints
    int32_t m,            // number of rows
    int32_t n             // number of columns
);

/* Purpose: Build a dense matrix from mpq_t input */
SLIP_info SLIP_build_dense_mpq
(
    SLIP_dense *A_output, // dense matrix, allocated but unused
    mpq_t **b,            // set of values as mpq_t
    int32_t m,            // number of rows
    int32_t n             // number of columns
);

/* Purpose: Build a dense matrix from mpfr_t input */
SLIP_info SLIP_build_dense_mpfr
(
    SLIP_dense *A_output, // Dense matrix, allocated but unused
    mpfr_t **b,           // Set of values as mpfr_t
    int32_t m,            // number of rows
    int32_t n,            // number of columns
    SLIP_options *option  // command options containing the prec for mpfr
);

//------------------------------------------------------------------------------
// double_mat: 2D double matrix
//------------------------------------------------------------------------------

// Creates a 2D double matrix, where A[i][j] is the (i,j)th entry.
// A[i] is a pointer to a row, of size n.

/* Purpose: This function creates a double matrix of size m*n. */
double** SLIP_create_double_mat
(
    int32_t m,     // number of rows
    int32_t n      // number of columns
);

/* Purpose: This function frees a dense double matrix.
 *
 * Input is a double*** mat and is destroyed on function completion.
 */
void SLIP_delete_double_mat
(
    double*** A,   // dense matrix
    int32_t m,     // number of rows of A
    int32_t n      // number of columns of A
);

//------------------------------------------------------------------------------
// int_mat: 2D int32 matrix
//------------------------------------------------------------------------------

// Creates a 2D int32 matrix, where A[i][j] is the (i,j)th entry.
// A[i] is a pointer to row i, of size n.

/* Purpose: This function creates an int matrix of size m*n. */
int32_t** SLIP_create_int_mat
(
    int32_t m,     // number of rows
    int32_t n      // number of columns
);

/* Purpose: This function deletes a dense int matrix.
 *
 * Input is a int*** mat and its dimensions.
 * Input mat is destroyed on function completion
 */
void SLIP_delete_int_mat
(
    int32_t*** A,  // dense matrix
    int32_t m,     // number of rows
    int32_t n      // number of columns
);

//------------------------------------------------------------------------------
// mpfr_mat:  a 2D mpfr_t matrix
//------------------------------------------------------------------------------

// Creates a 2D mpfr_t matrix, where A[i][j] is the (i,j)th entry.
// A[i] is a pointer to a row i, of size n.

/* Purpose: This function creates a mpfr_t matrix of size m*n with
 * precision prec
 */
mpfr_t** SLIP_create_mpfr_mat
(
    int32_t m,     // number of rows
    int32_t n,     // number of columns
    SLIP_options *option  // command options containing the prec for mpfr
);

/* Purpose: This function deletes a dense mpfr matrix.
 *
 * Input is a mpfr*** mat which is destroyed on completion
 */
void SLIP_delete_mpfr_mat
(
    mpfr_t ***A,   // Dense mpfr matrix
    int32_t m,     // number of rows of A
    int32_t n      // number of columns of A
);

//------------------------------------------------------------------------------
// mpq_mat:  a 2D mpq_t matrix
//------------------------------------------------------------------------------

// Creates a 2D mpq_t matrix, where A[i][j] is the (i,j)th entry.
// A[i] is a pointer to row i, of size n.

/* Purpose: This function creates a mpq_t matrix of size m*n. */
mpq_t** SLIP_create_mpq_mat
(
    int32_t m,     // number of rows
    int32_t n      // number of columns
);

/* Purpose: This function deletes a dense mpq matrix
 *
 * Input is a mpq_t*** matrix which is destroyed upon function completion
 */
void SLIP_delete_mpq_mat
(
    mpq_t***A,     // dense mpq matrix
    int32_t m,     // number of rows of A
    int32_t n      // number of columns of A
);

//------------------------------------------------------------------------------
// mpz_mat:  a 2D mpz_t matrix
//------------------------------------------------------------------------------

// Creates a 2D mpz_t matrix, where A[i][j] is the (i,j)th entry.
// A[i] is a pointer to row i, of size n.

/* Purpose: This function creates a dense mpz_t matrix of size m*n to
 * default size
 */
mpz_t** SLIP_create_mpz_mat
(
    int32_t m,     // number of rows
    int32_t n      // number of columns
);

/* Purpose: This function deletes a dense mpz matrix
 *
 * Input is a mpz_t*** matrix which is destoyed on completion
 */
void SLIP_delete_mpz_mat
(
    mpz_t ***A,     // The dense mpz matrix
    int32_t m,      // number of rows of A
    int32_t n       // number of columns of A
);

//------------------------------------------------------------------------------
// mpfr_vector: a 1D mpfr_t array
//------------------------------------------------------------------------------

// Creates a simple 1D array, where A[i] is an entry of type mpfr_t.

/* Purpose: This function creates a MPFR array of desired precision*/
mpfr_t* SLIP_create_mpfr_array
(
    int32_t n,     // size of the array
    SLIP_options *option  // command options containing the prec for mpfr
);

/* Purpose: This function clears the memory used for an mpfr array of size n.
 *
 * Input is a mpfr** array and its size. The input array is destroyed on output
 */
void SLIP_delete_mpfr_array
(
    mpfr_t** x,    // mpfr array to be deleted
    int32_t n      // size of x
);

//------------------------------------------------------------------------------
// mpq_vector: a 1D mpq_t array
//------------------------------------------------------------------------------

// Creates a simple 1D array, where A[i] is an entry of type mpq_t.

/* Purpose: This function creates an mpq array of size n.
 * This function must be called for all mpq arrays created.
 */
mpq_t* SLIP_create_mpq_array
(
    int32_t n      // size of the array
);

/* Purpose: This function clears the memory used for an mpq vector of size n.
 * Call this for all mpq vectors when done.
 *
 * Input is a mpq_t** array which is destroyed upon function completion
 */
void SLIP_delete_mpq_array
(
    mpq_t** x,     // mpq array to be deleted
    int32_t n      // size of x
);

//------------------------------------------------------------------------------
// mpz_vector: a 1D mpz_t array
//------------------------------------------------------------------------------

// Creates a simple 1D array, where A[i] is an entry of type mpz_t.

/* Purpose: This function creates an mpz array of size n and allocates
 * default size.
 */
mpz_t* SLIP_create_mpz_array
(
    int32_t n      // Size of x
);

/* Purpose: This function clears the memory used for an mpz vector of size n.
 * Call this function for all mpz vectors when done.
 *
 * Input is a mpz_t** array which is destroyed upon function completion
 *
 */
void SLIP_delete_mpz_array
(
    mpz_t ** x,     // mpz array to be deleted
    int32_t n       // Size of x
);

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// SLIP LU memory environment routines
//------------------------------------------------------------------------------

/*
 * Purpose: This function initializes the working evironment for SLIP LU
 * library.
 */
void SLIP_initialize (void);

/* Purpose: Initialize SLIP LU with user defined memory functions 
 */
void SLIP_initialize_expert
(
    void* (*MyMalloc) (size_t),                     // User defined malloc function
    void* (*MyRealloc) (void *, size_t, size_t),    // User defined realloc function
    void (*MyFree) (void*, size_t)                  // User defined free function
);

/*
 * Purpose: This function finalizes the working evironment for SLIP LU library.
 */
void SLIP_finalize (void);


//------------------------------------------------------------------------------
// Primary factorization & solve routines
//------------------------------------------------------------------------------

/*
 * Purpose: This function performs the symbolic ordering for SLIP LU. Currently,
 * there are four options: user defined order, COLAMD, AMD.
 */
SLIP_info SLIP_LU_analyze
(
    SLIP_LU_analysis *S,  // symbolic analysis (column permutation and nnz L,U)
    SLIP_sparse *A,       // Input matrix
    SLIP_options *option  // Control parameters
);

/* Purpose: This function performs the SLIP LU factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAQ = LDU
 * The determinant can be obtained as rhos[n-1]
 *
 *  A: input only, not modified
 *  L: allocated on input, modified on output
 *  U: allocated on input, modified on output
 *  S: input only, not modified
 *  rhos: allocated on input, modified on output
 *  pinv: allocated on input, modified on output
 *  option: input only, not modified
 */
SLIP_info SLIP_LU_factorize
(
    SLIP_sparse *L,         // lower triangular matrix
    SLIP_sparse *U,         // upper triangular matrix
    SLIP_sparse *A,         // matrix to be factored
    SLIP_LU_analysis *S,    // prior symbolic analysis
    mpz_t *rhos,            // sequence of pivots
    int32_t *pinv,          // inverse row permutation
    SLIP_options *option    // command options
);

SLIP_info SLIP_determinant
(
    SLIP_sparse *A,         // matrix that was factorized by SLIP_LU_factorize
    mpz_t *rhos,            // sequence of pivots from SLIP_LU_factorize
    mpq_t determinant       // determinant of A
) ;

// Solves Ax=b, returning the solution x as a double matrix
SLIP_info SLIP_solve_double
(
    double **x_doub,        // Solution vector stored as an double
    SLIP_sparse *A,         // Compressed column form full precision matrix A
    SLIP_LU_analysis *S,    // Column ordering
    SLIP_dense *b,          // Right hand side vectrors
    SLIP_options *option    // Control parameters
);

// Solves Ax=b, returning the solution x as an mpfr_t matrix
SLIP_info SLIP_solve_mpfr
(
    mpfr_t **x_mpfr,        // Solution vector stored as an mpfr_t array
    SLIP_sparse *A,         // Compressed column form full precision matrix A
    SLIP_LU_analysis *S,    // Column ordering
    SLIP_dense *b,          // Right hand side vectrors
    SLIP_options *option    // Control parameters
);

// Solves Ax=b, returning the solution x as an mpq_t matrix
SLIP_info SLIP_solve_mpq
(
    mpq_t **x_mpq,          // Solution vector stored as an mpq_t array
    SLIP_sparse *A,         // Compressed column form full precision matrix A
    SLIP_LU_analysis *S,    // Column ordering
    SLIP_dense *b,          // Right hand side vectrors
    SLIP_options *option    // Control parameters
);

/*
 * Purpose: This function permutes x to get it back in its original form.
 * That is x = Q*x.
 */
SLIP_info SLIP_permute_x
(
    mpq_t **x,            // Solution vector
    int32_t n,            // Size of solution vector
    int32_t numRHS,       // number of RHS vectors
    SLIP_LU_analysis *S   // symbolic analysis with the column ordering Q
);


/* Purpose: This function scales the x matrix if necessary */
SLIP_info SLIP_scale_x
(
    mpq_t **x,              // Solution matrix
    SLIP_sparse *A,         // matrix A
    SLIP_dense *b           // right hand side
);

/* Purpose: This function solves the linear system LD^(-1)U x = b.*/
SLIP_info SLIP_LU_solve     //solves the linear system LD^(-1)U x = b
(
    mpq_t **x,              // rational solution to the system
    SLIP_dense *b,          // right hand side vector
    mpz_t *rhos,            // sequence of pivots
    SLIP_sparse *L,         // lower triangular matrix
    SLIP_sparse *U,         // upper triangular matrix
    int32_t *pinv           // row permutation
);

// check and print a SLIP_sparse matrix
SLIP_info SLIP_spok  // returns a SLIP_LU status code
(
    SLIP_sparse *A,     // matrix to check
    int32_t print_level // 0: print nothing, 1: just errors, 2: terse, 3: all
) ;

/* Purpose: Convert the output mpq_t** solution vector obtained from
 * SLIP_Solve and SLIP_Permute_x from mpq_t** to double
 * x_doub has to be initialized before passed in
 */
SLIP_info SLIP_get_double_soln
(
    double **x_doub,      // double soln of size n*numRHS to Ax = b
    mpq_t  **x_mpq,       // mpq solution to Ax = b. x is of size n*numRHS
    int32_t n,            // Dimension of A, number of rows of x
    int32_t numRHS        // Number of right hand side vectors
) ;

/* Purpose: Convert the output mpq_t** solution vector obtained from
 * SLIP_Solve and SLIP_Permute_x from mpq_t** to mpfr_t**
 * x_mpfr has to be initialized before passed in
 */

SLIP_info SLIP_get_mpfr_soln
(
    mpfr_t **x_mpfr,      // mpfr solution of size n*numRHS to Ax = b
    mpq_t  **x_mpq,       // mpq solution of size n*numRHS to Ax = b.
    int32_t n,            // Dimension of A, number of rows of x
    int32_t numRHS        // Number of right hand side vectors
);

/*
 * Check the solution of the linear system
 * Performs a quick rational arithmetic A*x=b
 */
SLIP_info SLIP_check_solution
(
    SLIP_sparse *A,           // input matrix
    mpq_t **x,                // solution vector
    SLIP_dense *b             // right hand side
);

#endif
