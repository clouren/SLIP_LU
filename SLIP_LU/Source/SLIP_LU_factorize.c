//------------------------------------------------------------------------------
// SLIP_LU/SLIP_LU_factorize: exact sparse LU factorization
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs the SLIP LU factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAQ = LDU
 * The determinant of A can be obtained as determinant = rhos[n-1]
 *
 *  A: input only, not modified
 *  L: allocated on input, modified on output
 *  U: allocated on input, modified on output
 *  S: input only, not modified
 *  rhos: allocated on input, modified on output
 *  pinv: allocated on input, modified on output
 *  option: input only, not modified
 */

#define SLIP_FREE_WORKSPACE         \
    SLIP_delete_mpz_array(&x,n);    \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(col_loc);             \
    SLIP_FREE(pivs);                \
    SLIP_FREE(row_perm);            \
    SLIP_MPFR_CLEAR(temp);          \
    SLIP_MPZ_CLEAR(sigma);

# include "SLIP_LU_internal.h"

SLIP_info SLIP_LU_factorize
(
    SLIP_sparse *L,         // lower triangular matrix
    SLIP_sparse *U,         // upper triangular matrix
    SLIP_sparse *A,         // matrix to be factored
    SLIP_LU_analysis *S,    // stores guess on nnz and column permutation
    mpz_t * rhos,           // sequence of pivots
    int32_t* pinv,          // inverse row permutation
    SLIP_options* option    // command options
)
{
    // Input Check
    if (!A || !L || !U || !S || !rhos || !pinv || !option ||
        !A->p || !A->x || !A->i)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------
    // Begin timing factorization
    SLIP_info ok = SLIP_OK;
    int32_t n = A->n, k = 0, top, i, j, col, loc,
        lnz = 0, unz = 0, pivot, jnew, *xi = NULL, *h = NULL, *col_loc = NULL,
        *pivs = NULL, *row_perm = NULL;
    uint64_t size;
    mpz_t sigma; SLIP_MPZ_SET_NULL(sigma);
    mpfr_t temp; SLIP_MPFR_SET_NULL(temp);
    mpz_t* x = NULL ;

    SLIP_CHECK(SLIP_mpz_init(sigma));
    SLIP_CHECK(SLIP_mpfr_init2(temp, 256));
    // Sequence of chosen pivots
    pivs = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
    // Location of a column WRT the order
    col_loc = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
    // History vector
    h = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
    // Nonzero pattern
    xi = (int32_t*) SLIP_malloc(2*n* sizeof(int32_t));
    // Row permutation, inverse of pinv
    row_perm = (int32_t*) SLIP_malloc(n* sizeof(int32_t));

    if (!pivs || !col_loc || !h || !xi || !row_perm)
    {
        // out of memory: free everything and return
        SLIP_FREE_WORKSPACE ;
        return SLIP_OUT_OF_MEMORY;
    }
    slip_reset_int_array(pivs,n);
    slip_reset_int_array(h,n);

    //--------------------------------------------------------------------------
    // Compute a bound for the size of each entry in the matrix. This bound is
    // used to allocate the size of each entry in the x vector in order to
    // reduce the number of intermediate reallocations performed in the
    // triangular solve.
    //
    // This bound is based on a relaxation of Hadamard's bound
    //
    //--------------------------------------------------------------------------
    // Initialize sigma = largest entry in A

    SLIP_CHECK(SLIP_mpz_set(sigma, A->x[0]));

    // Get sigma = max(A)
    for (i = 1; i < A->nz; i++)
    {
        if(mpz_cmpabs(sigma,A->x[i]) < 0)
        {
            SLIP_CHECK(SLIP_mpz_set(sigma,A->x[i]));
        }
    }
    // sigma = |sigma|
    SLIP_CHECK(SLIP_mpz_abs(sigma,sigma));

    int32_t gamma = A->p[1];
    // get gamma as most dense column
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
    // Bound = gamma*log2(sigma sqrt(gamma))
    //--------------------------------------------------------------------------
    // temp = sigma*sqrt(gamma)
    SLIP_CHECK(SLIP_mpfr_mul_d(temp, temp, (double) sqrt(gamma), option->SLIP_MPFR_ROUND));
    // temp = log2(temp)
    SLIP_CHECK(SLIP_mpfr_log2(temp, temp, option->SLIP_MPFR_ROUND));
    // inner2 = temp
    double inner2;
    SLIP_CHECK(SLIP_mpfr_get_d(&inner2, temp, option->SLIP_MPFR_ROUND));
    // Free cache from log2. Even though mpfr_free_cache is called in
    // SLIP_LU_final(), it has to be called here to prevent memory leak in
    // some rare situations.
    SLIP_mpfr_free_cache();
    // bound = gamma * inner2+1
    int32_t bound = ceil(gamma*(inner2+1));
    // Ensure bound is at least 64 bit
    if (bound < 64) {bound = 64;}


    //--------------------------------------------------------------------------
    // Declare memory for x, L, and U
    //--------------------------------------------------------------------------

    // Initialize x
    x = slip_create_mpz_array2(n,bound);
    if (!x)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    // Initialize location based vectors
    for (i = 0; i < n; i++)
    {
        col_loc[S->q[i]] = i;
        pinv[i] = i;
        row_perm[i] = i;
    }

    // Allocate L and U.
    SLIP_CHECK(slip_sparse_alloc2(L, n, n, S->lnz));
    SLIP_CHECK(slip_sparse_alloc2(U, n, n, S->unz));

    //--------------------------------------------------------------------------
    // Iteration 0, must select pivot
    //--------------------------------------------------------------------------
    col = S->q[0];
    // x = A(:,col)
    SLIP_CHECK(slip_get_column(x, A, col)); 
    // top: nnz in column col
    top = n - ( (A->p[col+1]) - (A->p[col]) );
    j = 0;

    // Populate nonzero pattern
    for (i = A->p[col]; i < A->p[col+1]; i++)
    {
        xi[top+j] = A->i[i];
        j+=1;
    }
    // Get pivot
    SLIP_CHECK(slip_get_pivot(&pivot, x, pivs, n, top, xi, option->pivot,
        col, k, rhos, pinv, row_perm, option->tol));
    
    // Populate L and U
    for (j = top; j < n; j++)
    {
        jnew = xi[j];
        loc = pinv[jnew];

        //----------------------------------------------------------------------
        // U entries
        //----------------------------------------------------------------------
        if (loc <= k && ok == SLIP_OK)
        {
            // ith value of x[j]
            U->i[unz] = jnew;
            // Allocate memory for x[j]
            SLIP_CHECK(SLIP_mpz_sizeinbase(&size, x[jnew], 2));
            // GMP manual: Allocated size should be size+2
            SLIP_CHECK(SLIP_mpz_init2(U->x[unz],size+2));
            // Set U[x]
            SLIP_CHECK(SLIP_mpz_set(U->x[unz],x[jnew]));
            // Increment nnz of U
            unz++;
        }

        //----------------------------------------------------------------------
        // L entries
        //----------------------------------------------------------------------
        if (loc >= k && ok == SLIP_OK)
        {
            // ith value of x[j]
            L->i[lnz] = jnew;
            SLIP_CHECK(SLIP_mpz_sizeinbase(&size, x[jnew], 2));
            // GMP manual: Allocated size should be size+2
            SLIP_CHECK(SLIP_mpz_init2(L->x[lnz], size+2));
            // Set L[x]
            SLIP_CHECK(SLIP_mpz_set(L->x[lnz], x[jnew]));
            lnz++;
        }
    }

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
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
        SLIP_CHECK(slip_REF_triangular_solve(&top, L, A, k, xi, S->q, rhos,
            pinv, row_perm, col_loc, h, x));
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
            if (loc <= k && ok == SLIP_OK)
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
            if (loc >= k && ok == SLIP_OK)
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

    SLIP_FREE_WORKSPACE ;
    
    // Finalize L and U
    L->nz = lnz;
    U->nz = unz;
    // Finalize L->p, U->p
    L->p[n] = lnz;
    U->p[n] = unz;

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    if (ok == SLIP_OK)
    {
        // This cannot fail since the size of L and U are shrinking
        // Collapse L
        slip_sparse_collapse(L); 
        // Collapse U
        slip_sparse_collapse(U);
    }

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

    return ok;
}

