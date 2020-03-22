
// SLIP_type: an enum defined one of 5 types:
typedef enum
{
    SLIP_MPZ = 0 ;
    SLIP_MPQ = 1 ;
    SLIP_MPFR = 2 ;
    SLIP_INT64 = 3 ;
    SLIP_DOUBLE = 4 ;
}
SLIP_type ;

// SLIP_kind: an enum defined one of 5 types:
typedef enum
{
    SLIP_CSC = 0 ;
    SLIP_TRIPLET = 1 ;
    SLIP_DENSE = 2 ;
}
SLIP_kind ;

// SLIP_matrix:  a matrix of any type, in CSC, triplet, or dense format:
// 5 types * 3 kinds = 15 matrices.

typedef struct
{
    int64_t m ;     // Number of rows
    int64_t n ;     // Number of columns
    int64_t nzmax ; // Allocated size of A->i, A->j, and A->x
    int64_t nz ;    // Number of nonzeros in the matrix
    bool shallow ;  // if true, matrix is shallow.  When a shallow matrix
                    // is freed, the p, i, j, and x arrays are not modified.

    SLIP_kind kind ;    // CSC, triplet, or dense
    SLIP_type type ;    // mpz, mpq, mpfr, int64, or double

    int64_t *p;     // if CSC: Column pointers. Array size is n+1
                    // if triplet or dense: NULL.

    int64_t *i;     // if CSC or triplet: Row indices.
                    // Array size is nzmax, # of entries = nz.
                    // if dense: NULL

    int64_t *j;     // if triplet: Col indices. Array size is nzmax,
                    // # of entries = nz>
                    // if CSC or dense: j is NULL.

    // Values in matrix with size of nzmax, # of entries = nz:
    // Either a union, or a single (void *) pointer.
    union
    {
        mpz_t *mpz ;            // A->x.mpz
        mpq_t *mpq ;            // A->x.mpq
        mpfr_t *mpfr ;          // A->x.mpfr
        int64_t *int64 ;        // A->x.int64
        double *fp64 ;          // A->x.fp64
    } x ;

    mpq_t scale ; // scale factor for the matrix

} SLIP_matrix ;

// SLIP_matrix_allocate:  allocates an m-by-n matrix.
// if shallow is false: All components (p,j,i,x)  are allocated and set to zero
// if shallow is true:  All components (p,j,i,x) are NULL.
SLIP_Info SLIP_matrix_allocate
(
    SLIP_matrix **A,        // matrix to allocate
    int64_t m,              // # of rows
    int64_t n,              // # of columns
    int64_t nzmax,          // max # of entries
    SLIP_kind kind,         // CSC, triplet, or dense
    SLIP_type type,         // mpz, mpq, mpfr, int64, or double
    bool shallow,           // if true, matrix is shallow.
    SLIP_options *options
) ;

// SLIP_matrix_free
SLIP_Info SLIP_matrix_free
(
    SLIP_matrix **A,        // matrix to free
    SLIP_options *options
) ;

// SLIP_matrix_copy: makes a copy of a matrix (C not shallow, A might be)

    // 15x15
    // 3x3 functions (CSC, triplet, dense) <-> (CSC, triplet, dense)
    //  all 9 do any [5x5 typecasts]
    //      typecast function:  2 pointers, 2 SLIP_type, all 5x5, # entries

//  C->scale ?
//  A (double) to C (in mpz).  C->scale = ... ?
//  A (mpz) to C (in double).  use A->scale

SLIP_Info SLIP_matrix_copy
(
    SLIP_matrix **C,        // matrix to create (never shallow)
    SLIP_kind kind,         // CSC, triplet, or dense
    SLIP_type type,         // mpz_t, mpq_t, mpfr_t, int64_t, or double
    SLIP_Matrix *A,         // matrix to make a copy of (may be shallow)
    SLIP_options *option
) ;

// convert x in mpz (came from SLIP_solve ...)
SLIP_matrix *my_x, *x_soln ;
SLIP_solve (x_soln, ...) ;
SLIP_matrix_copy (&my_x, SLIP_SPARSE, SLIP_MPQ, x_soln, option) ;

// To access the value of the kth entry in any SLIP_matrix:
    SLIP_ENTRY (A, p, fp64) = x ;
    A->x.fp64 [p] = x

// To access a single entry in a dense SLIP_matrix:
#define INDEX(i,j,m) (i + j*A->m)

#define SLIPX(A,i,j,type)  A->x.type [i + j*(A->m)]
#define SLIPS(A,k,type)    A->x.type [k]

SLIPX (A, i, j, fp64)
SLIPXD (A, i, j)

#define SLIPXD(A,i,j)  SLIPX (A, i, j, fp64)
#define SLIPXZ(A,i,j)  SLIPX (A, i, j, mpz)
#define SLIPXQ(A,i,j)  SLIPX (A, i, j, mpq)
#define SLIPXR(A,i,j)  SLIPX (A, i, j, mpfr)

#define SLIPSD(A,k)  SLIPS (A, k, fp64)
#define SLIPSZ(A,k)  SLIPS (A, k, mpz)
#define SLIPSQ(A,k)  SLIPS (A, k, mpq)
#define SLIPSR(A,k)  SLIPS (A, k, mpfr)

    // A(i,j) = x
    SLIP_ENTRY (A, i, j, fp64) = x ;

    mpf (A->x [i][j]  ... )

    mpfr_set (A->x [i][j],  ... )
    mpfr_set (SLIPX (A,i,j,mpfr),  ... )

    mpfr_set (A->x [p],  ... )
    mpfr_set (A->x.mpfr [p],  ... )
    mpfr_set (SLIPS (A,p,mpfr),  ... )

    // A(i,j) = x
    SLIP (A, i, j, fp64) = x ;
    A->x.fp64 [i+j*(A->m)] = x

    // x = A(i,j) 
    x = SLIP (A, i, j, fp64) ;

_Generic ( )


