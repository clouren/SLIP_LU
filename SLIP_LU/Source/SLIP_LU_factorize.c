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
    SLIP_delete_mpz_array(&x,n);    \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(pivs);                \
    SLIP_FREE(row_perm);            \
    SLIP_MPFR_CLEAR(temp);          \
    SLIP_MPZ_CLEAR(sigma);

#define SLIP_FREE_WORKSPACE         \
    SLIP_FREE_WORK                  \
    SLIP_delete_sparse(&L);         \
    SLIP_delete_sparse(&U);         \
    SLIP_delete_mpz_array(&rhos,n); \
    SLIP_free(pinv);

#include "SLIP_LU_internal.h"

SLIP_info SLIP_LU_factorize
(
    // output:
    SLIP_sparse **L_handle, // lower triangular matrix
    SLIP_sparse **U_handle, // upper triangular matrix
    mpz_t **rhos_handle,    // sequence of pivots
    int32_t **pinv_handle,  // inverse row permutation
    // input:
    SLIP_sparse *A,         // matrix to be factored
    SLIP_LU_analysis *S,    // stores guess on nnz and column permutation
    SLIP_options* option    // command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_sparse *L = NULL ;
    SLIP_sparse *U = NULL ;
    mpz_t *rhos = NULL ;
    int32_t *pinv = NULL ;
    int32_t *xi = NULL ;
    int32_t *h = NULL ;
    int32_t *pivs = NULL ;
    int32_t *row_perm = NULL ;
    mpz_t* x = NULL ;

    mpz_t sigma; SLIP_MPZ_SET_NULL(sigma);
    mpfr_t temp; SLIP_MPFR_SET_NULL(temp);

    if (!L_handle || !U_handle || !rhos_handle || !pinv_handle)
    {
        return SLIP_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*U_handle) = NULL ;
    (*rhos_handle) = NULL ;
    (*pinv_handle) = NULL ;

    if (!A || !S || !option || !A->p || !A->x || !A->i)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    SLIP_info ok = SLIP_OK;
    int32_t n = A->n ;

    L = SLIP_create_sparse ( ) ;
    U = SLIP_create_sparse ( ) ;
    pinv = (int32_t *) SLIP_malloc (n * sizeof (int32_t)) ;
    rhos = SLIP_create_mpz_array (n) ;
    if (!L || !U || !pinv || !rhos)
    {
        // out of memory: free everything and return
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    int32_t k = 0, top, i, j, col, loc, lnz = 0, unz = 0, pivot, jnew ;
    size_t size ;

    SLIP_CHECK(SLIP_mpz_init(sigma));
    SLIP_CHECK(SLIP_mpfr_init2(temp, 256));

    // Indicator of which rows have been pivotal
    // pivs[i] = 1 if row i has been selected as a pivot
    // row, otherwise, pivs[i] < 0
    pivs = (int32_t*) SLIP_malloc(n* sizeof(int32_t));

    // h is the history vector utilized for the sparse REF
    // triangular solve algorithm. h serves as a global
    // vector which is repeatedly passed into the triangular
    // solve algorithm
    h = (int32_t*) SLIP_malloc(n* sizeof(int32_t));

    // xi is the global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L and U
    // for the triangular solve.
    xi = (int32_t*) SLIP_malloc(2*n* sizeof(int32_t));

    // Actual row permutation, the inverse of pinv. This
    // is used for sorting
    row_perm = (int32_t*) SLIP_malloc(n* sizeof(int32_t));

    if (!pivs || !h || !xi || !row_perm)
    {
        // out of memory: free everything and return
        SLIP_FREE_WORKSPACE ;
        return SLIP_OUT_OF_MEMORY;
    }
    // Initialize pivs and h; that is set pivs[i] = -1 and
    // h[i] = -1 for all i
    slip_reset_int_array(pivs,n);
    slip_reset_int_array(h,n);

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
    SLIP_CHECK(SLIP_mpz_set(sigma, A->x[0]));

    // Iterate throughout A and set sigma = max (A)
    for (i = 1; i < A->nz; i++)
    {
        if(mpz_cmpabs(sigma,A->x[i]) < 0)
        {
            SLIP_CHECK(SLIP_mpz_set(sigma,A->x[i]));
        }
    }

    // sigma = |sigma|
    SLIP_CHECK(SLIP_mpz_abs(sigma,sigma));

    // gamma is the number of nonzeros in the most dense column of A. First, we
    // initialize gamma to be the number of nonzeros in A(:,1).
    int32_t gamma = A->p[1];

    // Iterate throughout A and obtain gamma as the most dense column of A
    for (i = 1; i<n; i++)
    {
        if( gamma < A->p[i+1] - A->p[i])
        {
            gamma = A->p[i+1]-A->p[i];
        }
    }

    // temp = sigma
    SLIP_CHECK(SLIP_mpfr_set_z(temp, sigma, option->SLIP_MPFR_ROUND));

    //--------------------------------------------------------------------------
    // The bound is given as: gamma*log2(sigma sqrt(gamma))
    //--------------------------------------------------------------------------
    // temp = sigma*sqrt(gamma)
    SLIP_CHECK(SLIP_mpfr_mul_d(temp, temp, (double) sqrt(gamma),
        option->SLIP_MPFR_ROUND));
    // temp = log2(temp)
    SLIP_CHECK(SLIP_mpfr_log2(temp, temp, option->SLIP_MPFR_ROUND));
    // inner2 = temp
    double inner2;
    SLIP_CHECK(SLIP_mpfr_get_d(&inner2, temp, option->SLIP_MPFR_ROUND));
    // Free cache from log2. Even though mpfr_free_cache is called in
    // SLIP_LU_final(), it has to be called here to prevent memory leak in
    // some rare situations.
    SLIP_mpfr_free_cache();
    // bound = gamma * inner2+1. We add 1 to inner2 because log2(1) = 0
    int32_t bound = ceil(gamma*(inner2+1));
    // Ensure bound is at least 64 bit. In some rare cases the bound is very
    // small so we default to 64 bits.
    if (bound < 64) {bound = 64;}

    //--------------------------------------------------------------------------
    // Declare memory for x, L, and U
    //--------------------------------------------------------------------------

    // Initialize x. x is the global mpz_t vector of size n which is used
    // repeatedly in the sparse REF triangular solve.
    x = slip_create_mpz_array2(n,bound);
    if (!x)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    // Initialize location based vectors
    for (i = 0; i < n; i++)
    {
        pinv[i] = i;
        row_perm[i] = i;
    }

    // Allocate L and U.
    SLIP_CHECK(slip_sparse_alloc2(L, n, n, S->lnz));
    SLIP_CHECK(slip_sparse_alloc2(U, n, n, S->unz));

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
        //----------------------------------------------------------------------
        if (lnz + n > L->nzmax)
        {
            // Set L->nz = lnz
            L->nz = lnz;
            SLIP_CHECK(slip_sparse_realloc(L));
        }
        if (unz + n > U->nzmax)
        {
            // Set U->nz = unz
            U->nz = unz;
            SLIP_CHECK(slip_sparse_realloc(U));
        }

        // LDx = A(:,k)
        SLIP_CHECK(slip_REF_triangular_solve(&top, L, A, k, xi,
            (const int32_t *) (S->q),
            (const mpz_t   *) rhos,
            (const int32_t *) pinv,
            (const int32_t *) row_perm,
            h, x)) ;

        // Obtain pivot index
        SLIP_CHECK(slip_get_pivot(&pivot, x, pivs, n, top, xi, option->pivot,
            col, k, rhos, pinv, row_perm, option->tol));

        //----------------------------------------------------------------------
        // Iterate accross the nonzeros in x
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
                // Place the i location of the U->nz nonzero
                U->i[unz] = jnew;
                SLIP_CHECK(SLIP_mpz_sizeinbase(&size, x[jnew], 2));
                // GMP manual: Allocated size should be size+2
                SLIP_CHECK(SLIP_mpz_init2(U->x[unz], size+2));
                // Place the x value of the U->nz nonzero
                SLIP_CHECK(SLIP_mpz_set(U->x[unz], x[jnew]));
                // Increment U->nz
                unz++;
            }

            //------------------------------------------------------------------
            // loc >= k are rows below k, thus go to L
            //------------------------------------------------------------------
            if (loc >= k)
            {
                // Place the i location of the L->nz nonzero
                L->i[lnz] = jnew;
                SLIP_CHECK(SLIP_mpz_sizeinbase(&size, x[jnew], 2));
                // GMP manual: Allocated size should be size+2
                SLIP_CHECK(SLIP_mpz_init2(L->x[lnz], size+2));
                // Place the x value of the L->nz nonzero
                SLIP_CHECK(SLIP_mpz_set(L->x[lnz], x[jnew]));
                // Increment L->nz
                lnz++;
            }
        }
    }

    // Finalize L and U
    L->nz = lnz;
    U->nz = unz;
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
    for (i = 0; i < L->nz; i++)
    {
        L->i[i] = pinv[L->i[i]];
    }
    // Permute entries in U
    for (i = 0; i < U->nz; i++)
    {
        U->i[i] = pinv[U->i[i]];
    }

    //--------------------------------------------------------------------------
    // check the LU factorization (debugging only)
    //--------------------------------------------------------------------------

    #if 0
    SLIP_CHECK (SLIP_spok (L, option)) ;
    SLIP_CHECK (SLIP_spok (U, option)) ;
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

