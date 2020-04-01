//------------------------------------------------------------------------------
// SLIP_LU/SLIP_LU_factorize: exact sparse LU factorization
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs the SLIP LU factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAQ = LDU
 * The determinant of A can be obtained as determinant = rhos[n-1]
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

#define SLIP_FREE_WORK              \
    SLIP_matrix_free(&x, NULL);     \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(pivs);                \
    SLIP_FREE(row_perm);            \
    SLIP_MPFR_CLEAR(temp);          \
    SLIP_MPZ_CLEAR(sigma);

#define SLIP_FREE_ALL               \
    SLIP_FREE_WORK                  \
    SLIP_matrix_free(&L, NULL);     \
    SLIP_matrix_free(&U, NULL);     \
    SLIP_matrix_free(&rhos, NULL);  \
    SLIP_FREE(pinv);

#include "slip_LU_internal.h"

SLIP_info SLIP_LU_factorize
(
    // output:
    SLIP_matrix **L_handle,    // lower triangular matrix
    SLIP_matrix **U_handle,    // upper triangular matrix
    SLIP_matrix **rhos_handle, // sequence of pivots
    int64_t **pinv_handle,     // inverse row permutation
    // input:
    const SLIP_matrix *A,      // matrix to be factored
    const SLIP_LU_analysis *S, // stores guess on nnz and column permutation
    const SLIP_options* option // command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_REQUIRE (A, SLIP_CSC, SLIP_MPZ) ;

    if (!L_handle || !U_handle || !rhos_handle || !pinv_handle)
    {
        return SLIP_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*U_handle) = NULL ;
    (*rhos_handle) = NULL ;
    (*pinv_handle) = NULL ;

    if (!S  || !A->p || !A->i) 
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    SLIP_matrix *L = NULL ;    
    SLIP_matrix *U = NULL ;    
    SLIP_matrix *rhos = NULL ;       
    int64_t *pinv = NULL ;
    int64_t *xi = NULL ;
    int64_t *h = NULL ;
    int64_t *pivs = NULL ;
    int64_t *row_perm = NULL ;
    SLIP_matrix *x = NULL ;

    mpz_t sigma; SLIP_MPZ_SET_NULL(sigma);
    mpfr_t temp; SLIP_MPFR_SET_NULL(temp);

    SLIP_info info ;
    int64_t n = A->n ;
    pinv = (int64_t *) SLIP_malloc (n * sizeof (int64_t)) ;
    
    if (!pinv)
    {
        // out of memory: free everything and return
        SLIP_FREE_ALL ;
        return SLIP_OUT_OF_MEMORY;
    }

    int64_t k = 0, top, i, j, col, loc, lnz = 0, unz = 0, anz, pivot, jnew ;
    size_t size ;
    anz = A->p[n];
    
    SLIP_CHECK(SLIP_mpz_init(sigma));
    SLIP_CHECK(SLIP_mpfr_init2(temp, 256));

    // Indicator of which rows have been pivotal
    // pivs[i] = 1 if row i has been selected as a pivot
    // row, otherwise, pivs[i] < 0
    pivs = (int64_t*) SLIP_malloc(n* sizeof(int64_t));

    // h is the history vector utilized for the sparse REF
    // triangular solve algorithm. h serves as a global
    // vector which is repeatedly passed into the triangular
    // solve algorithm
    h = (int64_t*) SLIP_malloc(n* sizeof(int64_t));

    // xi is the global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L and U
    // for the triangular solve.
    xi = (int64_t*) SLIP_malloc(2*n* sizeof(int64_t));

    // Actual row permutation, the inverse of pinv. This
    // is used for sorting
    row_perm = (int64_t*) SLIP_malloc(n* sizeof(int64_t));

    if (!pivs || !h || !xi || !row_perm)
    {
        // out of memory: free everything and return
        SLIP_FREE_ALL  ;
        return SLIP_OUT_OF_MEMORY;
    }
    // Initialize pivs and h; that is set pivs[i] = -1 and
    // h[i] = -1 for all i
    for (i = 0; i < n; i++)
    {
        h[i] = -1;
        pivs[i] = -1;
    }

    //--------------------------------------------------------------------------
    // This section of the code computes a bound for the worst case bit-length
    // of each entry in the matrix. This bound is used to allocate the size of
    // each number in the global x vector. As a result of this allocation,
    // computing the values in L and U via repeated triangular solves will
    // not require intermediate memory reallocations from the GMP library.
    //
    // This bound is based on a relaxation of sparse Hadamard's bound
    //--------------------------------------------------------------------------

    // sigma is the largest initial entry in A. First we initialize sigma to be
    // the first nonzero in A
    SLIP_CHECK(SLIP_mpz_set(sigma, A->x.mpz[0]));

    // Iterate throughout A and set sigma = max (A)
    for (i = 1; i < anz; i++)
    {
        int r ;
        SLIP_CHECK (SLIP_mpz_cmpabs (&r, sigma, A->x.mpz [i])) ;
        if(r < 0)
        {
            SLIP_CHECK(SLIP_mpz_set(sigma,A->x.mpz[i]));
        }
    }

    // sigma = |sigma|
    SLIP_CHECK(SLIP_mpz_abs(sigma,sigma));

    // gamma is the number of nonzeros in the most dense column of A. First, we
    // initialize gamma to be the number of nonzeros in A(:,1).
    int64_t gamma = A->p[1];

    // Iterate throughout A and obtain gamma as the most dense column of A
    for (i = 1; i<n; i++)
    {
        if( gamma < A->p[i+1] - A->p[i])
        {
            gamma = A->p[i+1]-A->p[i];
        }
    }

    // temp = sigma
    SLIP_CHECK(SLIP_mpfr_set_z(temp, sigma, SLIP_GET_ROUND(option)));

    //--------------------------------------------------------------------------
    // The bound is given as: gamma*log2(sigma sqrt(gamma))
    //--------------------------------------------------------------------------

    // temp = sigma*sqrt(gamma)
    SLIP_CHECK(SLIP_mpfr_mul_d(temp, temp, (double)sqrt(gamma), SLIP_GET_ROUND(option)));
    // temp = log2(temp)
    SLIP_CHECK(SLIP_mpfr_log2(temp, temp, SLIP_GET_ROUND(option)));
    // inner2 = temp
    double inner2;
    SLIP_CHECK(SLIP_mpfr_get_d(&inner2, temp, SLIP_GET_ROUND(option)));
    // Free cache from log2. Even though mpfr_free_cache is called in
    // SLIP_LU_final(), it has to be called here to prevent memory leak in
    // some rare situations due to the usage of the mpfr_log2 function.
    // as per MPFR's documentation
    SLIP_mpfr_free_cache();
    // bound = gamma * inner2+1. We add 1 to inner2 because log2(1) = 0
    int64_t bound = ceil(gamma*(inner2+1));
    // Ensure bound is at least 64 bit. In some rare cases the bound is very
    // small so we default to 64 bits.
    if (bound < 64) {bound = 64;}

    //--------------------------------------------------------------------------
    // Declare memory for x, rhos, L, and U
    //--------------------------------------------------------------------------

    // Create rhos, a global dense mpz_t matrix of dimension n*1
    SLIP_CHECK (SLIP_matrix_allocate(&rhos, SLIP_DENSE, SLIP_MPZ, n, 1, n,
        false, true, option));

    // Create x, a global dense mpz_t matrix of dimension n*1. Unlike rhos, the
    // second boolean parameter, x->init is set to false to avoid initializing
    // each mpz entry of x with default size, which should be bound calculated
    // above as gamma*log2(sigma sqrt(gamma)) instead
    SLIP_CHECK (SLIP_matrix_allocate(&x, SLIP_DENSE, SLIP_MPZ, n, 1, n,
        false, false, option));

    // Allocate L and U without initializing each entry.
    // L and U are allocated to have nnz(L) = the estimate from symbolic 
    // analysis. However, unlike traditional matrix allocation, the second
    // boolean parameter here is set to false, so the individual values of 
    // L and U are not allocated. Instead, a more efficient method to 
    // allocate these values is done in the factorization to reduce 
    // memory usage.
    SLIP_CHECK (SLIP_matrix_allocate(&L, SLIP_CSC, SLIP_MPZ, n, n, S->lnz,
        false, false, option));
    SLIP_CHECK (SLIP_matrix_allocate(&U, SLIP_CSC, SLIP_MPZ, n, n, S->unz,
        false, false, option));

    for (i = 0; i < n; i++)
    {
        // Allocate memory for entries of x
        SLIP_CHECK(SLIP_mpz_init2(x->x.mpz[i], bound));

        // Initialize location based vectors
        pinv[i] = i;
        row_perm[i] = i;
    }

    //--------------------------------------------------------------------------
    // Iterations 0:n-1 (1:n in standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // Column pointers for column k of L and U
        L->p[k] = lnz;
        U->p[k] = unz;
        col = S->q[k];

        //----------------------------------------------------------------------
        // Reallocate memory if necessary
        // if lnz+n > L->nzmax, L needs to expand to accomodate new nonzeros.
        // To do so, we double the size of the L and U matrices.
        //----------------------------------------------------------------------
        if (lnz + n > L->nzmax)
        {
            // Double the size of L
            SLIP_CHECK(slip_sparse_realloc(L));
        }
        if (unz + n > U->nzmax)
        {
            // Double the size of U
            SLIP_CHECK(slip_sparse_realloc(U));
        }

        //----------------------------------------------------------------------
        // Triangular solve to compute LDx = A(:,k)
        //----------------------------------------------------------------------
        SLIP_CHECK(slip_ref_triangular_solve(&top, L, A, k, xi,
            (const int64_t *) (S->q),
            rhos,
            (const int64_t *) pinv,
            (const int64_t *) row_perm,
            h, x)) ;

        //----------------------------------------------------------------------
        // Obtain pivot
        //----------------------------------------------------------------------
        SLIP_CHECK(slip_get_pivot(&pivot, x, pivs, n, top, xi, option->pivot,
            col, k, rhos, pinv, row_perm, option->tol));

        //----------------------------------------------------------------------
        // Populate L and U. We iterate across all nonzeros in x
        //----------------------------------------------------------------------
        for (j = top; j < n; j++)
        {
            jnew = xi[j];
            // Location of x[j] in final matrix
            loc = pinv[jnew];

            //------------------------------------------------------------------
            // loc <= k are rows above k, thus go to U
            //------------------------------------------------------------------
            if (loc <= k)
            {
                // Place the i location of the unz nonzero
                U->i[unz] = jnew;
                // Find the size in bits of x[j]
                SLIP_CHECK(SLIP_mpz_sizeinbase(&size, x->x.mpz[jnew], 2));
                // GMP manual: Allocated size should be size+2
                SLIP_CHECK(SLIP_mpz_init2(U->x.mpz[unz], size+2));
                // Place the x value of the unz nonzero
                SLIP_CHECK(SLIP_mpz_set(U->x.mpz[unz], x->x.mpz[jnew]));
                // Increment unz
                unz++;
            }

            //------------------------------------------------------------------
            // loc >= k are rows below k, thus go to L
            //------------------------------------------------------------------
            if (loc >= k)
            {
                // Place the i location of the lnz nonzero
                L->i[lnz] = jnew;
                // Set the size of x[j]
                SLIP_CHECK(SLIP_mpz_sizeinbase(&size, x->x.mpz[jnew], 2));
                // GMP manual: Allocated size should be size+2
                SLIP_CHECK(SLIP_mpz_init2(L->x.mpz[lnz], size+2));
                // Place the x value of the lnz nonzero
                SLIP_CHECK(SLIP_mpz_set(L->x.mpz[lnz], x->x.mpz[jnew]));
                // Increment lnz
                lnz++;
            }
        }
    }

    // Finalize L->p, U->p
    L->p[n] = lnz;
    U->p[n] = unz;
    
    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    // free everything, but keep L, U, rhos, and pinv
    SLIP_FREE_WORK ;

    // This cannot fail since the size of L and U are shrinking.
    // Collapse L
    slip_sparse_collapse(L);
    // Collapse U
    slip_sparse_collapse(U);

    //--------------------------------------------------------------------------
    // finalize the row indices in L and U
    //--------------------------------------------------------------------------

    // Permute entries in L
    for (i = 0; i < lnz; i++)
    {
        L->i[i] = pinv[L->i[i]];
    }
    // Permute entries in U
    for (i = 0; i < unz; i++)
    {
        U->i[i] = pinv[U->i[i]];
    }

    //--------------------------------------------------------------------------
    // check the LU factorization (debugging only)
    //--------------------------------------------------------------------------

    #if 0
    SLIP_CHECK (SLIP_matrix_check (L, option)) ;
    SLIP_CHECK (SLIP_matrix_check (U, option)) ;
    #endif

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    (*L_handle) = L ;
    (*U_handle) = U ;
    (*rhos_handle) = rhos ;
    (*pinv_handle) = pinv ;
    return (SLIP_OK) ;
}

