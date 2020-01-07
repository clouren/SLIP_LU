//------------------------------------------------------------------------------
// SLIP_LU/SLIP_LU_internal: include file for internal use in SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SLIP_LU.h instead.

#ifndef SLIP_LU_internal
#define SLIP_LU_internal

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------C Libraries------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Standard C libraries
# include <setjmp.h>
# include <stdbool.h>
# include <assert.h>
# include <stdlib.h>
# include <math.h>
# include <stdio.h>
# include <stdint.h>
# include <string.h>
# include <time.h>
# include <stdarg.h>

// SuiteSparse headers
# include "SuiteSparse_config.h"
# include "colamd.h"
# include "amd.h"

// gmp libraries
# include <gmp.h>
# include <mpfr.h>

#include "SLIP_LU.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Common Macros----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


#ifdef MATLAB_MEX_FILE

    # include "mex.h"
    # include "matrix.h"
    // use the MATLAB memory manager
    #define SLIP_MEMORY_MALLOC  mxMalloc
    #define SLIP_MEMORY_CALLOC  mxCalloc
    #define SLIP_MEMORY_REALLOC mxRealloc
    #define SLIP_MEMORY_FREE    mxFree

#else

    // use the ANSI C memory manager
    #define SLIP_MEMORY_MALLOC  malloc
    #define SLIP_MEMORY_CALLOC  calloc
    #define SLIP_MEMORY_REALLOC realloc
    #define SLIP_MEMORY_FREE    free

#endif


#define SLIP_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SLIP_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SLIP_FLIP(i) (-(i)-2)
#define SLIP_UNFLIP(i) (((i) < 0) ? SLIP_FLIP(i) : (i))
#define SLIP_MARKED(Ap,j) (Ap [j] < 0)
#define SLIP_MARK(Ap,j) { Ap [j] = SLIP_FLIP (Ap [j]) ; }

// SLIP_CHECK(method) is a macro that calls a slip method and checks the
// status; if a failure occurs, it prints the detailed error message, frees all
// allocated workspace and returns the error status to the caller.
// To use SLIP_CHECK, the #include'ing file must declare a SLIP_info ok,
// and must define SLIP_FREE_WORKSPACE as a macro that frees all workspace
// if an error occurs. The method can be a scalar ok as well, so that
// SLIP_CHECK(ok) works.

// the default is to free nothing
#ifndef SLIP_FREE_WORKSPACE
#define SLIP_FREE_WORKSPACE
#endif

#define SLIP_CHECK(method)                                              \
{                                                                       \
    ok = method ;                                                       \
    if (ok != SLIP_OK)                                                  \
    {                                                                   \
        SLIP_FREE_WORKSPACE ;                                           \
        return (ok) ;                                                   \
    }                                                                   \
}

#ifdef SLIP_LU_TCOV
    /* include this header to redefine SLIP_MEMORY_REALLOC (used in SLIP_gmp.c)
     * for memory test and to use macro GOTCHA */
    #include "../Tcov/tcov_malloc_test.h"
#endif

// Tolerance used in the pivoting schemes. This number can be anything in
// between 0 and 1. A value of 0 selects the diagonal element exclusively and a
// value of 1 selects the smallest or largest pivot exclusively only in a
// tolerance pivoting based method
#define SLIP_DEFAULT_TOL 0.1

// Check parameter. If this = 1 then the solution to the system is checked
// for accuracy
#define SLIP_DEFAULT_CHECK false

// Pivoting scheme used for SLIP LU.
//  SLIP_SMALLEST = 0,              Smallest pivot
//  SLIP_DIAGONAL = 1,              Diagonal pivoting
//  SLIP_FIRST_NONZERO = 2,         First nonzero per column chosen as pivot
//  SLIP_TOL_SMALLEST = 3,          Diagonal pivoting with tolerance for small
//  SLIP_TOL_LARGEST = 4,           Diagonal pivoting with tolerance for large
//  SLIP_LARGEST = 5                Largest pivot
#define SLIP_DEFAULT_PIVOT SLIP_TOL_SMALLEST

// Column ordering used.
//  SLIP_NO_ORDERING = 0,           None: Not recommended for sparse matrices
//  SLIP_COLAMD = 1,                COLAMD: Default
//  SLIP_AMD = 2                    AMD
#define SLIP_DEFAULT_ORDER SLIP_COLAMD

// Defines printing to be done
#define SLIP_DEFAULT_PRINT_LEVEL 0

// MPFR precision used (quad is default)
#define SLIP_DEFAULT_PRECISION 128

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
// 
// SLIP LU utilizes MPFR_RNDN by default. 

#define SLIP_DEFAULT_MPFR_ROUND MPFR_RNDN;

// Size of mpz_t, mpq_t and mpfr_t values
#define SIZE_MPZ  sizeof(mpz_t)
#define SIZE_MPQ  sizeof(mpq_t)
#define SIZE_MPFR sizeof(mpfr_t)


// Field access macros for MPZ/MPQ/MPFR struct
// (similar definition in gmp-impl.h and mpfr-impl.h)

#define MPZ_SIZ(x)   ((x)->_mp_size)
#define MPZ_PTR(x)   ((x)->_mp_d)
#define MPZ_ALLOC(x) ((x)->_mp_alloc)
#define MPQ_NUM(x)   mpq_numref(x)
#define MPQ_DEN(x)   mpq_denref(x)
#define MPFR_MANT(x) ((x)->_mpfr_d)
#define MPFR_EXP(x)  ((x)->_mpfr_exp)
#define MPFR_PREC(x) ((x)->_mpfr_prec)
#define MPFR_SIGN(x) ((x)->_mpfr_sign)
#define MPFR_REAL_PTR(x) (&((x)->_mpfr_d[-1])) /*re-define but same result*/
/* Invalid exponent value (to track bugs...) */
#define MPFR_EXP_INVALID \
 ((mpfr_exp_t) 1 << (GMP_NUMB_BITS*sizeof(mpfr_exp_t)/sizeof(mp_limb_t)-2))

/* Macros to set the pointer in mpz_t/mpq_t/mpfr_t variable to NULL. It is best
 * practice to call these macros immediately after mpz_t/mpq_t/mpfr_t variable
 * is declared, and before the mp*_init function is called. It would help to
 * prevent error when SLIP_MP*_CLEAR is called before the variable is
 * successfully initialized.
 */

#define SLIP_MPZ_SET_NULL(x)                \
    MPZ_PTR(x) = NULL;                      \
    MPZ_SIZ(x) = 0;                         \
    MPZ_ALLOC(x) = 0;

#define SLIP_MPQ_SET_NULL(x)                \
    MPZ_PTR(MPQ_NUM(x)) = NULL;             \
    MPZ_SIZ(MPQ_NUM(x)) = 0;                \
    MPZ_ALLOC(MPQ_NUM(x)) = 0;              \
    MPZ_PTR(MPQ_DEN(x)) = NULL;             \
    MPZ_SIZ(MPQ_DEN(x)) = 0;                \
    MPZ_ALLOC(MPQ_DEN(x)) = 0;

#define SLIP_MPFR_SET_NULL(x)               \
    MPFR_MANT(x) = NULL;                    \
    MPFR_PREC(x) = 0;                       \
    MPFR_SIGN(x) = 1;                       \
    MPFR_EXP(x) = MPFR_EXP_INVALID;

/* GMP does not give a mechanism to tell a user when an mpz, mpq, or mpfr
 * item has been cleared; thus, if mp*_clear is called on an object that
 * has already been cleared, gmp will crash. It is also not possible to
 * set a mp*_t = NULL. Thus, this mechanism modifies the internal GMP
 * size of entries to avoid crashing in the case that a mp*_t is cleared
 * multiple times.
 */

#define SLIP_MPZ_CLEAR(x)                   \
{                                           \
    if ((x) != NULL && MPZ_PTR(x) != NULL)  \
    {                                       \
        mpz_clear(x);                       \
        SLIP_MPZ_SET_NULL(x);               \
    }                                       \
}

#define SLIP_MPQ_CLEAR(x)                   \
{                                           \
    SLIP_MPZ_CLEAR(MPQ_NUM(x));             \
    SLIP_MPZ_CLEAR(MPQ_DEN(x));             \
}

#define SLIP_MPFR_CLEAR(x)                  \
{                                           \
    if ((x) != NULL && MPFR_MANT(x) != NULL)\
    {                                       \
        mpfr_clear(x);                      \
        SLIP_MPFR_SET_NULL(x);              \
    }                                       \
}


// ============================================================================
//                           Internal Functions
// ============================================================================

/* Purpose: This function creates an mpz array of size n and allocates
 * memory for numbers of bit size prec. If the relative size of numbers is
 * known ahead of time, this is more efficient than the
 * SLIP_create_mpz_array
 */

mpz_t* slip_create_mpz_array2
(
    int32_t n,     // size of the array
    int32_t size   // Relative size of numbers
);

/*
 * Purpose: This function resets an mpz array of size n with the nonzero pattern
 * given. This is more efficient than iterating accross all nonzeros in vector x
 */
SLIP_info slip_reset_mpz_array
(
    mpz_t* x,      // mpz array to be reset
    int32_t n,     // size of x
    int32_t top,   // beginning of the nonzero pattern
    int32_t* xi    // nonzero pattern
);


/*
 * Purpose: This function initializes an int vector of size n and sets the value
 * equal to -1. This function is used for the history and pivot vectors.
 */
SLIP_info slip_reset_int_array
(
    int32_t* h,           // int vector to be reset
    int32_t n             // size of the int vector
);

/*
 * Purpose: This function resets an int vector of size n and sets each term
 * equal to -1 with nonzero pattern given. This is more efficient than resetting
 * each term individually
 */
SLIP_info slip_reset_int_array2
(
    int32_t* h,    // int vector to be reset
    int32_t n,     // size of h
    int32_t top,   // beginning of nonzero pattern
    int32_t* xi    // nonzero pattern
);

/* Purpose: This function takes as input a mpz_t** array and divides it by a
 * mpz_t constant storing the solution in a mpq_t** array. This is used
 * internally to divide the solution vector by the determinant of the matrix.
 *
 * On output, the contents of the array x2 are modified
 *
 */
SLIP_info slip_array_div // divides the x vector by the determinant
(
    mpq_t** x2,     // solution of x/det
    mpz_t** x,      // input vector
    mpz_t det,      // given determinant of matrix
    int32_t n,      // size of x and x2
    int32_t numRHS  // number of rhs vectors
);

/*
 * Purpose: This function multiplies vector x by the determinant of matrix.
 *
 * On output the contents of the x vector is modified
 */
SLIP_info slip_array_mul // multiplies vector x by the determinant of matrix
(
    mpz_t** x,      // matrix to be multiplied
    mpz_t det,      // given determinant of matrix
    int32_t n,      // size of x
    int32_t numRHS  // number of RHS vectors
);

/* Purpose: This function performs sparse REF forward substitution This is
 * essentially the same as the sparse REF triangular solve applied to each
 * column of the right hand side vectors. Like the normal one, this
 * function expects that the vector x is dense. As a result,the nonzero
 * pattern is not computed and each nonzero in x is iterated across.
 * The system to solve is LDx = x
 *
 * Out output, the mpz_t** x structure is modified
 *
 */
SLIP_info slip_forward_sub
(
    SLIP_sparse *L,   // lower triangular matrix
    mpz_t **x,        // right hand side matrix of size n*numRHS
    mpz_t *rhos,      // sequence of pivots used in factorization
    int32_t numRHS    // number of columns in x
);

/* Purpose: This function performs sparse REF backward substitution. In essense
 * it solves the sysem Ux = x. Note that prior to this, we expect x to be
 * multiplied by the determinant of A.
 *
 * The input argument x is modified on output
 *
 */

SLIP_info slip_back_sub  // performs sparse REF backward substitution
(
    SLIP_sparse *U,   // input upper triangular matrix
    mpz_t **bx,       // right hand side matrix of size n*numRHS
    int32_t numRHS    // number of columns in bx
)  ;

/*
 * Purpose: p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] in
 * to c
 * This function is lightly modified from CSparse
 *
 */
SLIP_info slip_cumsum
(
    int32_t *p,      // vector to store the sum of c
    int32_t *c,      // vector which is summed
    int32_t n        // size of c
);

/* Purpose: This function performs a depth first search of the graph of the
 * matrix starting at node j. The output of this function is the set of nonzero
 * indices in the xi vector
 *
 */
void slip_dfs // performs a dfs of the graph of the matrix starting at node j
(
    int32_t *top,    // beginning of stack
    int32_t j,       // What node to start DFS at
    SLIP_sparse* L,  // matrix which represents the Graph of L
    int32_t* xi,     // the nonzero pattern
    int32_t* pstack, // workspace vector
    int32_t* pinv    // row permutation
);

/* Purpose: This function converts a double array of size n to an appropriate
 * mpz array of size n. To do this, the number is multiplied by 10^17 then, the
 * GCD is found. This function allows the use of matrices in double precision to
 * work with SLIP LU
 * NOTE: First element of input double array must be nonzero
 */
SLIP_info slip_expand_double_array
(
    mpz_t *x_out,//integral final array
    double* x,  //double array that needs to be made integral
    mpq_t scale,//the scaling factor used (x_out = scale * x)
    int32_t n,   //size of x
    SLIP_options* option
);

/* Purpose: This function converts a double matrix of size m*n to an appropriate
 * mpz array of size m*n. To do this, the number is multiplied by 10^17 then,
 * the GCD is found. This function allows the use of matrices in double
 * precision to work with SLIP LU
 * NOTE: First element of input double matrix must be nonzero
 */
SLIP_info slip_expand_double_mat
(
    mpz_t **x_out,   //integral final mat
    double** x,     // double matrix that needs to be made integral
    mpq_t scale,    // the scaling factor used (x_out = scale * x)
    int32_t m,      // row number of x
    int32_t n,       // column number of x
    SLIP_options* option
);

/* Purpose: This function converts a mpfr array of size n and precision prec to
 * an appropriate mpz array of size n. To do this, the number is multiplied by
 * the appropriate power of 10 then the gcd is found. This function allows mpfr
 * arrays to be used within SLIP LU
 * NOTE: First element of input mpfr_t array must be nonzero
 */
SLIP_info slip_expand_mpfr_array
(
    mpz_t *x_out,   //integral final array
    mpfr_t* x,      // mpfr array to be expanded
    mpq_t scale,    // scaling factor used (x_out = scale*x)
    int32_t n,      // size of x
    SLIP_options *option// command options containing the prec for mpfr
);

/* Purpose: This function converts a mpfr matrix of size m*n and precision prec
 * to an appropriate mpz matrix of size m*n. To do this, the number is
 * multiplied by the appropriate power of 10 then the gcd is found. This
 * function allows mpfr arrays to be used within SLIP LU
 * NOTE: First element of input mpfr_t matrix must be nonzero
 */
SLIP_info slip_expand_mpfr_mat
(
    mpz_t **x_out,     //integral final mat
    mpfr_t** x,        // mpfr matrix to be expanded
    mpq_t scale,       // scaling factor used (x_out = scale*x)
    int32_t m,         // size of x
    int32_t n,         // size of x
    SLIP_options *option// command options containing the prec for mpfr
);

/* Purpose: This function converts a mpq array of size n into an appropriate mpz
 * array of size n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq arrays to be used in SLIP LU
 */
SLIP_info slip_expand_mpq_array
(
    mpz_t *x_out, //integral final array
    mpq_t* x,     //mpq array that needs to be converted
    mpq_t scale,  //scaling factor. x_out = scale*x
    int32_t n     //size of x
);

/* Purpose: This function converts a mpq matrix of size m*n into an appropriate
 * mpz matrix of size m*n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq matrix to be used in SLIP LU
 *
 * on output, x2 is modified
 *
 */
SLIP_info slip_expand_mpq_mat
(
    mpz_t **x_out,//integral final mat
    mpq_t **x,    // mpq mat that needs to be converted
    mpq_t scale,  // scaling factor. x_out = scale*x
    int32_t m,    // number of rows of x
    int32_t n     // number of columns of x
);

/* Purpose: This function obtains column k from matrix A and stores it in the
 * dense vector x
 *
 * On exit, x either contains the kth column of A or is NULL
 */
SLIP_info slip_get_column //extract k-th column from A, i.e., x=A(:,k)
(
    mpz_t* x,       // A(:,k)
    SLIP_sparse* A, // input matrix
    int32_t k       // column to extract
) ;


/* This function performs the pivoting for the SLIP LU factorization.
 * The optional Order is:
 *     0: Smallest pivot 
 *     1: Natural/Diagonal pivoting
 *     2: Choose first nonzero
 *     3: Diagonal with tolerance and smallest pivot (default)
 *     4: Diagonal with tolerance and largest pivoting
 *     5: Largest pivot
 *
 * On output, the pivs, rhos, pinv, and row_perm arrays are all modified
 *
 */

SLIP_info slip_get_pivot
(
    int32_t *pivot,
    mpz_t* x,       // kth column of L and U
    int32_t* pivs,  // vecor indicating which rows have been pivotal
    int32_t n,      // dimension of the problem
    int32_t top,    // nonzero pattern is located in xi[top..n-1]
    int32_t* xi,    // nonzero pattern of x
    int32_t order,  // what kind of pivoting to use (see above description)
    int32_t col,    // current column of A (real kth column i.e., q[k])
    int32_t k,      // iteration of the algorithm
    mpz_t* rhos,    // vector of pivots
    int32_t* pinv,  // row permutation
    int32_t* row_perm,// opposite of pinv. if pinv[i] = j then row_perm[j] = i
    double tolerance// tolerance used if some tolerance based pivoting is used
) ;

/* Purpose: This function selects the pivot element as the largest in the column
 * This is activated if the user sets option->pivot = SLIP_LARGEST
 * NOTE: This pivoting scheme is NOT recommended for SLIP LU
 *
 * On output, the index of the largest pivot is returned
 *
 */
SLIP_info slip_get_largest_pivot
(
    int32_t *pivot,  // index of largest pivot
    mpz_t* x,        // kth column of L and U
    int32_t* pivs,   // vector which indicates whether each row has been pivotal
    int32_t n,       // dimension of problem
    int32_t top,     // nonzero pattern is located in xi[top..n-1]
    int32_t* xi      // nonzero pattern of x
);

/* This function obtains the first eligible nonzero pivot
 * This is enabled if the user sets option->pivot = SLIP_FIRST_NONZERO
 * NOTE: This pivoting scheme is not recommended
 *
 * On output, the kth pivot is returned
 *
 */
SLIP_info slip_get_nonzero_pivot // find the first eligible nonzero pivot
(
    int32_t *pivot,   // index of nonzero pivot
    mpz_t* x,         // kth column of L and U
    int32_t* pivs,    // vector indicating which rows are pivotal
    int32_t n,        // size of x
    int32_t top,      // nonzero pattern is located in xi[top..n-1]
    int32_t* xi       // nonzero pattern of x
);

/* Purpose: This function selects the pivot element as the smallest in the
 * column. This is activated by default or if the user sets
 * option->pivot = SLIP_SMALLEST
 * NOTE: This is the recommended pivoting scheme for SLIP LU
 *
 * On output, the index of kth pivot is returned
 *
 */
SLIP_info slip_get_smallest_pivot
(
    int32_t *pivot, // index of smallest pivot
    mpz_t* x,       // kth column of L and U
    int32_t* pivs,  // vector indicating whether each row has been pivotal
    int32_t n,      // dimension of problem
    int32_t top,    // nonzeros are stored in xi[top..n-1]
    int32_t* xi     // nonzero pattern of x
);

/* Purpose: This function print the basic info about SLIP_LU library*/
void slip_lu_info(void);

/*
 * Purpose: This function allocates a SLIP LU matrix of size n*m with array size
 * nzmax. This function manually initializes each entry in A->x therefore they
 * are immediately ready to be operated on. This is less efficient but more user
 * friendly.
 */
SLIP_info slip_sparse_alloc
(
    SLIP_sparse* A,// sparse matrix data structure to be allocated
    int32_t n,     // number of columns
    int32_t m,     // number of rows (recall m=n assumed)
    int32_t nzmax  // size of allocated i and x arrays
);

/*
 * Purpose: This function allocates a SLIP LU matrix of size n*m with array size
 * nzmax. This version does not allocate individual the values in x. As a
 * result, it is more memory efficient, but also less user friendly.
 */
SLIP_info slip_sparse_alloc2
(
    SLIP_sparse* A, // sparse matrix data structure to be allocated
    int32_t n,      // number of columns
    int32_t m,      // number of rows (recall m=n assumed)
    int32_t nzmax   // size of allocated i and x arrays
);

/*
 * Purpose: This function collapses a SLIP matrix. Essentially it shrinks the
 * size of x and i. so that they only take up the number of elements in the
 * matrix. For example if A->nzmax = 1000 but A->nz = 500, r and x are of size
 * 1000, so this function shrinks them to size 500
 */
SLIP_info slip_sparse_collapse
(
    SLIP_sparse* A // matrix to be shrunk
);

/*
 * Purpose: This function expands a SLIP LU matrix by doubling its size. It
 * merely expands x and i and does not initialize/allocate the values!
 */
SLIP_info slip_sparse_realloc
(
    SLIP_sparse* A // the matrix to be expanded
);

/*
 * Purpose: This function allocates a SLIP dense matrix of size n*m.
 * This function manually initializes each entry in A->x therefore they
 * are immediately ready to be operated on.
 */
SLIP_info slip_dense_alloc
(
    SLIP_dense* A, // dense matrix data structure to be allocated
    int32_t m,     // number of rows
    int32_t n      // number of columns
);

/*
 * Purpose: This function populates the SLIP_sparse A by the ccf
 * stored vectors i, p, and x
 */
SLIP_info SLIP_mpz_populate_mat
(
    SLIP_sparse* A,   // matrix to be populated
    int32_t* I,       // row indices
    int32_t* p,       // column pointers
    mpz_t* x,         // set of values in A
    int32_t n,        // size of the matrix A
    int32_t nz        // number of nonzeros in the matrix A
);

/*
 * Purpose: This function computes the reach of column k of A on the graph of L
 * mathematically that is: xi = Reach(A(:,k))_G_L
 */
void slip_reach    // compute the reach of column k of A on the graph of L
(
    int32_t *top,
    SLIP_sparse* L,   // matrix representing graph of L
    SLIP_sparse* A,   // input matrix
    int32_t k,        // column of A of interest
    int32_t* xi,      // nonzero pattern
    int32_t* pinv     // row permutation
)  ;

/*
 * Purpose: This function sorts the xi vector with respect to the current row
 * permutation. This sort is efficient as its complexity is |x| log |x|.
 * The idea of the sort is that you have xi[top, top+1, ...]. We essentially
 * mask them and sort the masked vector (which is with respect to the row
 * permutation). We then unmask them to get the correct value. For instance, the
 * correct sorted order could be [4 2 5 1] because of the column permutation.
 */

void slip_sort_xi
(
    int32_t* xi,        // nonzero pattern
    int32_t top,        // nonzeros are stored in xi[top..n-1]
    int32_t n,          // size of problem
    int32_t* pinv,      // inverse row permutation
    int32_t* row_perm   // opposite of pinv. if pinv[j] = k then row_perm[k] = j
);


/*
 * Purpose: This function performs the sparse REF triangular solve. i.e.,
 * (LD) x = A(:,k). The algorithm is described in the paper; however in essence
 * it computes the nonzero pattern xi, then performs a sequence of IPGE
 * operations on the nonzeros to obtain their final value. All operations are
 * gauranteed to be integral. There are various enhancements in this code used
 * to reduce the overall cost of the operations and minimize operations as much
 * as possible.
 */
SLIP_info slip_REF_triangular_solve // performs the sparse REF triangular solve
(
    int32_t *top_output,      // Output the beginning of nonzero pattern
    SLIP_sparse* L,           // partial L matrix
    SLIP_sparse* A,           // input matrix
    int32_t k,                // iteration of algorithm
    int32_t* xi,              // nonzero pattern vector
    int32_t* q,               // column permutation
    mpz_t* rhos,              // sequence of pivots
    int32_t* pinv,            // inverse row permutation
    int32_t* row_perm,        // row permutation
    int32_t* col_loc,         // column permutation
    int32_t* h,               // history vector
    mpz_t* x                  // solution of system ==> kth column of L and U
);

/*
 * Purpose: This function converts triplet matrix into compressed column
 * matrix A
 */
SLIP_info slip_trip_to_mat
(
    SLIP_sparse *A,     //matrix stored in ccf that will take B
    int32_t *I,         // Row indices.
    int32_t *J,         // Column indices
    mpz_t *x,           // Values in the matrix
    int32_t n,          // Dimension of the matrix
    int32_t nz          // Number of nonzeros in the matrix
);

#endif
