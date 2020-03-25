//------------------------------------------------------------------------------
// SLIP_LU/Include/SLIP_LU.h: user #include file for SLIP_LU.
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#ifndef SLIP_LU_H
#define SLIP_LU_H

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
//---------------------Include files required by SLIP LU------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>

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
// Error codes
//------------------------------------------------------------------------------

// Most SLIP_LU functions return a code that indicates if it was successful
// or not. Otherwise the code returns a pointer to the object that was created
// or it returns void (in the case that an object was deleted)

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
    SLIP_SMALLEST = 0,      // Smallest pivot
    SLIP_DIAGONAL = 1,      // Diagonal pivoting
    SLIP_FIRST_NONZERO = 2, // First nonzero per column chosen as pivot
    SLIP_TOL_SMALLEST = 3,  // Diagonal pivoting with tol for smallest pivot.
                            //   (Default)
    SLIP_TOL_LARGEST = 4,   // Diagonal pivoting with tol. for largest pivot
    SLIP_LARGEST = 5        // Largest pivot
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

// This struct serves as a global struct to define all options

typedef struct SLIP_options
{
    SLIP_pivot pivot;     // Type of pivoting scheme used.
    SLIP_col_order order; // Type of column ordering scheme used
    double tol;           // User specified tolerance for SLIP_TOL_SMALLEST and
                          // SLIP_TOL_LARGEST
    int print_level;      // 0: print nothing, 1: just errors,
                          // 2: terse, with basic stats from COLAMD/AMD and SLIP
                          // 3: all, with matrices and results
    uint64_t prec;        // Precision used to output file if MPFR is chosen
    mpfr_rnd_t round;     // Type of MPFR rounding used
    bool check;           // Set true if the solution to the system should be checked
} SLIP_options;

/* Purpose: Create and return SLIP_options pointer with default parameters
 * upon successful allocation, which are defined in SLIP_LU_internal.h
 * To free it, simply use SLIP_FREE(option)
 */
SLIP_options* SLIP_create_default_options(void);

////////////////////////////////////////////////////////////////////////////////
// new SLIP_matrix-based methods
////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------
// SLIP_matrix: a sparse CSC, sparse triplet, or dense matrix
//------------------------------------------------------------------------------

// SLIP LU uses a single matrix data type, SLIP_matrix, which can be held in
// one of three kinds of formats:  sparse CSC (compressed sparse column),
// sparse triplet, and dense:

typedef enum
{
    SLIP_CSC = 0,           // matrix is in compressed sparse column format
    SLIP_TRIPLET = 1,       // matrix is in sparse triplet format
    SLIP_DENSE = 2          // matrix is in dense format
}
SLIP_kind ;

// Each of the three formats can have values of 5 different data types: mpz_t,
// mpq_t, mpfr_t, int64_t, and double:

typedef enum
{
    SLIP_MPZ = 0,           // matrix of mpz_t integers
    SLIP_MPQ = 1,           // matrix of mpq_t rational numbers
    SLIP_MPFR = 2,          // matrix of mpfr_t
    SLIP_INT64 = 3,         // matrix of int64_t integers
    SLIP_FP64 = 4           // matrix of doubles
}
SLIP_type ;

// This gives a total of 15 different matrix types.  Not all functions accept
// all 15 matrices types, however.

// Suppose A is an m-by-n matrix with nz <= nzmax entries.  These terms are
// A->m A->n A->nz and A->nzmax.  The p, i, j, and x components are defined as:

// (0) SLIP_CSC:  A sparse matrix in CSC (compressed sparse column) format.
//      A->p is an int64_t array of size n+1, A->i is an int64_t array of size
//      nzmax (with nz <= nzmax), and A->x.type is an array of size nzmax of
//      matrix entries ('type' is one of mpz, mpq, mpfr, int64, or fp64).  The
//      row indices of column j appear in A->i [A->p [j] ... A->p [j+1]-1], and
//      the values appear in the same locations in A->x.type.  The A->j array
//      is NULL.

// (1) SLIP_TRIPLET:  A sparse matrix in triplet format.  A->i and A->j are
//      both int64_t arrays of size nzmax, and A->x.type is an array of values
//      of the same size.  The kth tuple has row index A->i [k], column index
//      A->j [k], and value A->x.type [k], with 0 <= k < nz.  The A->p array
//      is NULL.

// (2) SLIP_DENSE:  A dense matrix.  The integer arrays A->p, A->i, and A->j
//      are all NULL.  A->x.type is a pointer to an array of size m*n, stored
//      in column-oriented format.  The value of A(i,j) is A->x.type [p]
//      with p = i + j*A->m.

// The SLIP_matrix may contain 'shallow' components, A->p, A->i, A->j, and
// A->x.  For example, if A->p_shallow is true, then a non-NULL A->p is a
// pointer to a read-only array, and the A->p array is not freed by
// SLIP_matrix_free.  If A->p is NULL (for a triplet or dense matrix), then
// A->p_shallow has no effect.

typedef struct
{
    int64_t m ;         // number of rows
    int64_t n ;         // number of columns
    int64_t nzmax ;     // size of A->i, A->j, and A->x
    int64_t nz ;        // # nonzeros in the matrix 
                        // For tripet and CSC matrices,
                        // all user visible functions will 
                        // have A->nz = A->nzmax. For dense
                        // matrices, A->nz is never used.
    SLIP_kind kind ;    // CSC, triplet, or dense
    SLIP_type type ;    // mpz, mpq, mpfr, int64, or fp64 (double)

    int64_t *p ;        // if CSC: column pointers, an array size is n+1.
                        // if triplet or dense: A->p is NULL.
    bool p_shallow ;    // if true, A->p is shallow.

    int64_t *i ;        // if CSC or triplet: row indices, of size nzmax.
                        // if dense: A->i is NULL.
    bool i_shallow ;    // if true, A->i is shallow.

    int64_t *j ;        // if triplet: column indices, of size nzmax.
                        // if CSC or dense: A->j is NULL.
    bool j_shallow ;    // if true, A->j is shallow.

    union               // A->x.type has size nzmax.
    {
        mpz_t *mpz ;            // A->x.mpz
        mpq_t *mpq ;            // A->x.mpq
        mpfr_t *mpfr ;          // A->x.mpfr
        int64_t *int64 ;        // A->x.int64
        double *fp64 ;          // A->x.fp64
    } x ;
    bool x_shallow ;    // if true, A->x.type is shallow.

    mpq_t scale ;   // scale factor for mpz matrices (never shallow)
                        // For all matrices who's type is not mpz,
                        // mpz_scale = 1. 

} SLIP_matrix ;

//------------------------------------------------------------------------------
// SLIP_matrix_allocate: allocate an m-by-n SLIP_matrix
//------------------------------------------------------------------------------

// if shallow is false: All components (p,i,j,x) are allocated and set to zero,
//                      and then shallow flags are all false.

// if shallow is true:  All components (p,i,j,x) are NULL, and their shallow
//                      flags are all true.  The user can then set A->p,
//                      A->i, A->j, and/or A->x accordingly, from their own
//                      arrays.

SLIP_info SLIP_matrix_allocate
(
    SLIP_matrix **A_handle, // matrix to allocate
    SLIP_kind kind,         // CSC, triplet, or dense
    SLIP_type type,         // mpz, mpq, mpfr, int64, or double
    int64_t m,              // # of rows
    int64_t n,              // # of columns
    int64_t nzmax,          // max # of entries
    bool shallow,           // if true, matrix is shallow.  A->p, A->i, A->j,
                            // A->x are all returned as NULL and must be set
                            // by the caller.  All A->*_shallow are returned
                            // as true.
    bool init,              // If true, and the data types are mpz, mpq, or
                            // mpfr, the entries are initialized (using the 
                            // appropriate SLIP_mp*_init function). If false,
                            // the mpz, mpq, and mpfr arrays are malloced but not 
                            // initialized. Utilized internally to reduce memory
    SLIP_options *option
) ;

//------------------------------------------------------------------------------
// SLIP_matrix_free: free a SLIP_matrix
//------------------------------------------------------------------------------

SLIP_info SLIP_matrix_free
(
    SLIP_matrix **A_handle, // matrix to free
    SLIP_options *option
) ;

//------------------------------------------------------------------------------
// SLIP_matrix_nnz: # of entries in a matrix
//------------------------------------------------------------------------------

int64_t SLIP_matrix_nnz     // return # of entries in A, or -1 on error
(
    SLIP_matrix *A,         // matrix to query
    SLIP_options *option
) ;

//------------------------------------------------------------------------------
// SLIP_matrix_copy: makes a copy of a matrix
//------------------------------------------------------------------------------

// SLIP_matrix_copy: make a copy of a SLIP_matrix, into another kind and type.

SLIP_info SLIP_matrix_copy
(
    SLIP_matrix **C,        // matrix to create (never shallow)
    // inputs, not modified:
    SLIP_kind kind,         // CSC, triplet, or dense
    SLIP_type type,         // mpz_t, mpq_t, mpfr_t, int64_t, or double
    SLIP_matrix *A,         // matrix to make a copy of (may be shallow)
    SLIP_options *option
) ;

//------------------------------------------------------------------------------
// SLIP_matrix macros
//------------------------------------------------------------------------------

// These macros simplify the access to entries in a SLIP_matrix.
// The type parameter is one of: mpq, mpz, mpfr, int64, or fp64.

// To access the kth entry in a SLIP_matrix using 1D linear addressing,
// in any matrix kind (CSC, triplet, or dense), in any type:
#define SLIP_1D(A,k,type) ((A)->x.type [k])

// To access the (i,j)th entry in a 2D SLIP_matrix, in any type:
#define SLIP_2D(A,i,j,type) SLIP_1D (A, (i)+(j)*((A)->m), type)

//------------------------------------------------------------------------------
// SLIP_LU_analysis: symbolic pre-analysis
//------------------------------------------------------------------------------

/* This struct stores the column permutation for LU and the guess on nnz for
 * L and U */

typedef struct
{
    int64_t *q;     // Column permutation for LU
    int64_t lnz;    // Approximate number of nonzeros in L.
                    // i.e., initial size of L
    int64_t unz;    // Approximate number of nonzeros in U.
                    // i.e., initial size of U
} SLIP_LU_analysis;

// The symbolic analysis object is created by SLIP_LU_analyze.

/* Purpose: This function frees the memory of the SLIP_LU_analysis struct
 *
 * Input is the SLIP_LU_analysis structure, it is destroyed on function
 * termination.
 */
void SLIP_LU_analysis_free        
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
    void* (*MyMalloc) (size_t),                     // User defined malloc
    void* (*MyRealloc) (void *, size_t, size_t),    // User defined realloc
    void (*MyFree) (void*, size_t)                  // User defined free
);

/*
 * Purpose: This function finalizes the working evironment for SLIP LU library.
 */
void SLIP_finalize (void);


//------------------------------------------------------------------------------
// Primary factorization & solve routines
//------------------------------------------------------------------------------

/* Purpose: SLIP_LU_analyze performs the symbolic ordering for SLIP LU.
 * Currently, * there are three options: user defined order, COLAMD, AMD.
 */
SLIP_info SLIP_LU_analyze
(
    SLIP_LU_analysis **S, // symbolic analysis (column permutation and nnz L,U)
    SLIP_matrix *A,       // Input matrix
    SLIP_options *option  // Control parameters
);

/* Purpose: SLIP_LU_factorize performs the SLIP LU factorization. This
 * factorization is done via n iterations of the sparse REF triangular solve
 * function. The overall factorization is PAQ = LDU.  The determinant can be
 * obtained as rhos->x.mpz[n-1]
 *
 *  L: undefined on input, created on output
 *  U: undefined on input, created on output
 *  rhos: undefined on input, created on output
 *  pinv: undefined on input, created on output
 *
 *  A: input only, not modified
 *  S: input only, not modified
 *  option: input only, not modified
 */
SLIP_info SLIP_LU_factorize
(
    // output:
    SLIP_matrix **L_handle,     // lower triangular matrix
    SLIP_matrix **U_handle,     // upper triangular matrix
    SLIP_matrix **rhos_handle,  // sequence of pivots
    int64_t **pinv_handle,      // inverse row permutation
    // input:
    SLIP_matrix *A,             // matrix to be factored
    SLIP_LU_analysis *S,        // stores guess on nnz and column permutation
    SLIP_options* option        // command options
) ;


/* Purpose: Solve the linear system Ax = b. This is the simplest way to 
 * use the SLIP LU package. This function encompasses both factorization
 * and solve and returns the solution vector in the user desired type. 
 * It can be thought of as an exact version of MATLAB sparse backslash.
 */
SLIP_info SLIP_backslash
(
    // Output
    SLIP_matrix **X_handle, // Final solution vector
    // Input
    SLIP_type type,         // Type of output desired
                            // Must be SLIP_MPQ, SLIP_MPFR,
                            // or SLIP_FP64
    SLIP_matrix *A,         // Input matrix
    SLIP_matrix *b,         // Right hand side vector(s)
    SLIP_options* option    // Command options 
);

/* Purpose: SLIP_LU_solve solves the linear system LD^(-1)U x = b.*/
SLIP_info SLIP_LU_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SLIP_matrix **X_handle,  // rational solution to the system
    // input:
    SLIP_matrix *b,         // right hand side vector
    const SLIP_matrix *A,   // Input matrix
    const SLIP_matrix *L,   // lower triangular matrix
    const SLIP_matrix *U,   // upper triangular matrix
    const SLIP_matrix *rhos,// sequence of pivots
    SLIP_LU_analysis *S,    // symbolic analysis struct
    const int64_t *pinv,    // row permutation
    bool check,             // Set true to check solution (for debugging)
    SLIP_options* option    // Command options
);

// SLIP_matrix_check: check and print a SLIP_sparse matrix
SLIP_info SLIP_matrix_check  // returns a SLIP_LU status code
(
    SLIP_matrix *A,       // matrix to check
    SLIP_options* option  // Determine print level
) ;


//------------------------------------------------------------------------------
//---------------------------SLIP GMP/MPFR Functions----------------------------
//------------------------------------------------------------------------------

/* The following functions are the SLIP LU interface to the GMP/MPFR libary.
 * Each corresponding GMP/MPFR function is given a wrapper to ensure that no
 * memory leaks or crashes occur. All covered GMP functions can be found in
 * SLIP_gmp.c
 *
 */

// The GMP library does not handle out-of-memory failures.  However, it does
// provide a mechanism for passing function pointers that replace GMP's use of
// malloc, realloc, and free.  This mechanism is used to provide a try/catch
// mechanism for memory allocation errors, using setjmp and longjmp.

// When a GMP function is called, this wrapper keeps track of a list of objects
// allocated by that function.  The list is started fresh each time a GMP
// function is called.  If any allocation fails, the NULL pointer is not
// returned to GMP.  Instead, all allocated blocks in the list are freed,
// and slip_gmp_allocate returns directly to wrapper.

#ifdef SLIP_LU_TCOV
    /* include this header to redefine SLIP_MEMORY_REALLOC (used in SLIP_gmp.c)
     * for memory test and to use macro GOTCHA */
    #include "../Tcov/tcov_malloc_test.h"
#endif

SLIP_info SLIP_gmp_fprintf(FILE *fp, const char *format, ... );

SLIP_info SLIP_gmp_printf(const char *format, ... );

SLIP_info SLIP_gmp_fscanf(FILE *fp, const char *format, ... );

#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info SLIP_mpfr_asprintf (char **str, const char *template, ... );

SLIP_info SLIP_mpfr_free_str (char *str);
#endif

SLIP_info SLIP_mpfr_fprintf(FILE *fp, const char *format, ... );

SLIP_info SLIP_mpz_init(mpz_t x) ;

SLIP_info SLIP_mpz_init2(mpz_t x, const size_t size) ;

SLIP_info SLIP_mpz_set(mpz_t x, const mpz_t y);

SLIP_info SLIP_mpz_set_ui(mpz_t x, const uint64_t y) ;

SLIP_info SLIP_mpz_set_si(mpz_t x, const int64_t y);

SLIP_info SLIP_mpz_set_d(mpz_t x, const double y);

SLIP_info SLIP_mpz_get_d(double *x, const mpz_t y) ;

SLIP_info SLIP_mpz_get_si(int64_t *x, const mpz_t y);

SLIP_info SLIP_mpz_set_q(mpz_t x, const mpq_t y) ;

SLIP_info SLIP_mpz_mul(mpz_t a, const mpz_t b, const mpz_t c) ;

#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info SLIP_mpz_add(mpz_t a, const mpz_t b, const mpz_t c);

SLIP_info SLIP_mpz_addmul(mpz_t x, const mpz_t y, const mpz_t z) ;
#endif

SLIP_info SLIP_mpz_submul(mpz_t x, const mpz_t y, const mpz_t z);

SLIP_info SLIP_mpz_divexact(mpz_t x, const mpz_t y, const mpz_t z) ;

SLIP_info SLIP_mpz_gcd(mpz_t x, const mpz_t y, const mpz_t z) ;

SLIP_info SLIP_mpz_lcm(mpz_t lcm, const mpz_t x, const mpz_t y) ;

SLIP_info SLIP_mpz_abs(mpz_t x, const mpz_t y) ;

SLIP_info SLIP_mpz_cmp(int *r, const mpz_t x, const mpz_t y) ;

SLIP_info SLIP_mpz_cmpabs(int *r, const mpz_t x, const mpz_t y) ;

SLIP_info SLIP_mpz_cmp_ui(int *r, const mpz_t x, const uint64_t y) ;

SLIP_info SLIP_mpz_sgn(int *sgn, const mpz_t x) ;

SLIP_info SLIP_mpz_sizeinbase(size_t *size, const mpz_t x, int64_t base) ;

SLIP_info SLIP_mpq_init(mpq_t x) ;

SLIP_info SLIP_mpq_set(mpq_t x, const mpq_t y);

SLIP_info SLIP_mpq_set_z(mpq_t x, const mpz_t y);

SLIP_info SLIP_mpq_set_d(mpq_t x, const double y) ;

SLIP_info SLIP_mpq_set_ui(mpq_t x, const uint64_t y, const uint64_t z);

SLIP_info SLIP_mpq_set_si (mpq_t x, const int64_t y, const uint64_t z) ;

SLIP_info SLIP_mpq_set_num(mpq_t x, const mpz_t y) ;

SLIP_info SLIP_mpq_set_den(mpq_t x, const mpz_t y);

SLIP_info SLIP_mpq_get_den(mpz_t x, const mpq_t y) ;

SLIP_info SLIP_mpq_get_d(double *x, const mpq_t y) ;

SLIP_info SLIP_mpq_abs(mpq_t x, const mpq_t y) ;

SLIP_info SLIP_mpq_add(mpq_t x, const mpq_t y, const mpq_t z) ;

SLIP_info SLIP_mpq_mul(mpq_t x, const mpq_t y, const mpq_t z) ;

SLIP_info SLIP_mpq_div(mpq_t x, const mpq_t y, const mpq_t z) ;

SLIP_info SLIP_mpq_cmp(int *r, const mpq_t x, const mpq_t y) ;

SLIP_info SLIP_mpq_cmp_ui(int *r, const mpq_t x,
                    const uint64_t num, const uint64_t den) ;

SLIP_info SLIP_mpq_sgn(int *sgn, const mpq_t x) ;

SLIP_info SLIP_mpq_equal(int *r, const mpq_t x, const mpq_t y) ;

SLIP_info SLIP_mpfr_init2(mpfr_t x, const uint64_t size) ;

SLIP_info SLIP_mpfr_set(mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_set_d(mpfr_t x, const double y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_set_si(mpfr_t x, int64_t y, const mpfr_rnd_t rnd);

SLIP_info SLIP_mpfr_set_q(mpfr_t x, const mpq_t y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_set_z(mpfr_t x, const mpz_t y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_get_z(mpz_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_get_q(mpq_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_get_d(double *x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_get_si(int64_t *x, const mpfr_t y, const mpfr_rnd_t rnd);

SLIP_info SLIP_mpfr_mul(mpfr_t x, const mpfr_t y, const mpfr_t z,
                    const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_mul_d(mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_div_d(mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_ui_pow_ui(mpfr_t x, const uint64_t y, const uint64_t z,
                    const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_log2(mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SLIP_info SLIP_mpfr_sgn(int *sgn, const mpfr_t x) ;

SLIP_info SLIP_mpfr_free_cache(void);

#endif
