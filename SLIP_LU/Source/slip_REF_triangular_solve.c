//------------------------------------------------------------------------------
// SLIP_LU/slip_REF_triangular_solve: sparse REF triangular solve
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/*
 * Purpose: This function performs the sparse REF triangular solve. i.e.,
 * (LD) x = A(:,k). The algorithm is described in the paper; however in essence
 * it computes the nonzero pattern xi, then performs a sequence of IPGE
 * operations on the nonzeros to obtain their final value. All operations are
 * gauranteed to be integral. There are various enhancements in this code used
 * to reduce the overall cost of the operations and minimize operations as much
 * as possible.
 */

/*

    q       an array of size n+1, where j = q [k] if the original column j
            is the kth column of the LU factorization.  Unmodified.

    pinv

    row_perm

    xi      workspace of size 2n, uninitialized on input.  On output,
            xi [top..n-1] contain the nonzero pattern of L(:,k), in
            strictly sorted order.  The row indices in xi [top..n-1] reflect
            the row indices in the original A matrix, not the final L and U,
            since the row permutation.

    top_output  undefined on input.  On output, top = (*top_output) gives the
            size of the pattern xi [top..n-1].

    x       workspace of size n, uninitialized on input.  On output,
            x [i] is the value of L(i,k), where i is in the pattern
            xi [top..n-1].  Other entries of x are undefined on output.

    h       workspace of size n,
            uninitialized on input and undefined on output.

*/


#include "SLIP_LU_internal.h"

SLIP_info slip_REF_triangular_solve // performs the sparse REF triangular solve
(
    int32_t *top_output,      // Output the beginning of nonzero pattern
    SLIP_sparse* L,           // partial L matrix
    SLIP_sparse* A,           // input matrix
    int32_t k,                // constructing L(:,k)
    int32_t* xi,              // nonzero pattern vector
    const int32_t* q,         // column permutation, not modified
    mpz_t* rhos,              // sequence of pivots
    const int32_t* pinv,      // inverse row permutation
    const int32_t* row_perm,  // row permutation
    int32_t* h,               // history vector
    mpz_t* x                  // solution of system ==> kth column of L and U
)
{
    // inputs have been validated in SLIP_LU_factorize.c
    int32_t j, jnew, i, inew, p, m, n, col, sgn, top ;
    SLIP_info ok;

    //--------------------------------------------------------------------------
    // Initialize REF TS by getting nonzero patern of x && obtaining A(:,k)
    //--------------------------------------------------------------------------

    // Size of matrix and the dense vectors
    n = A->n;

    // This triangular solve is computing L(:,k) and U(:,k) for the
    // factorization L*D*U = A*Q.  The column permutation is held in the
    // permutation vector q, and col = q [k] is the column of A that has been
    // permuted to the kth column of A*Q.
    col = q[k];

    // Obtain nonzero pattern of L(:,k) in xi[top..n]
    slip_reach(&top, L, A, col, xi, pinv);

    // Sort xi [top..n-1] wrt sequence of pivots
    slip_sort_xi(xi, top, n, pinv, row_perm);

    // Reset x[i] = 0 for all i in nonzero pattern xi [top..n-1]
    SLIP_CHECK(slip_reset_mpz_array(x, n, top, xi));

    // Set x[col] = 0.  A(col,col) is the diagonal entry in the original
    // matrix.  The pivot search prefers to select the diagonal, if it is
    // present.  It may be not present in A itself, at all, and thus also
    // not in the pattern xi [top..n-1].  The value x[col] is set to zero
    // here, in case the entry A(col,col) is not present, so that the pivot
    // search query the value of the diagonal.
    SLIP_CHECK(slip_mpz_set_ui(x[col], 0));

    // Reset h[i] = -1 for all i in nonzero pattern xi [top..n-1]
    SLIP_CHECK(slip_reset_int_array2(h, n, top, xi));

    // Set x = A(:,q(k))
    SLIP_CHECK(slip_get_column(x, A, col));

    //--------------------------------------------------------------------------
    // Iterate across nonzeros in x
    //--------------------------------------------------------------------------
    for (p = top; p < n; p++)
    {
        /* Finalize x[j] */
        j = xi[p];                         // First nonzero term
        jnew = pinv[j];                    // Location of nonzero term
        // Check if x[j] == 0, if so continue to next nonzero
        SLIP_CHECK(SLIP_mpz_sgn(&sgn, x[j]));
        if (sgn == 0) {continue;}          // x[j] = 0 no work must be done

        // x[j] is nonzero
        if (jnew < k)                      // jnew < k implies entries in U
        {
            //------------------------------------------------------------------
            // History update: Bring x[j] to its final value
            //------------------------------------------------------------------
            if (h[j] < jnew - 1)           // HU must be performed
            {
                // x[j] = x[j] * rho[j-1]
                SLIP_CHECK(SLIP_mpz_mul(x[j],x[j],rhos[jnew-1]));

                if (h[j] > -1)
                {
                    // x[j] = x[j] / rho[h[j]]
                    SLIP_CHECK(SLIP_mpz_divexact(x[j],x[j],rhos[h[j]]));
                }
            }

            //------------------------------------------------------------------
            // IPGE updates
            //------------------------------------------------------------------

            // ----------- Iterate across nonzeros in Lij ---------------------
            for (m = L->p[jnew]; m < L->p[jnew+1]; m++)
            {
                i = L->i[m];               // i value of Lij
                inew = pinv[i];            // i location of Lij
                if (inew > jnew)
                {
                    /*************** If lij==0 then no update******************/
                    SLIP_CHECK(SLIP_mpz_sgn(&sgn, L->x[m]));
                    if (sgn == 0) {continue;}

                    // lij is nonzero. Check if x[i] is nonzero
                    SLIP_CHECK(SLIP_mpz_sgn(&sgn, x[i]));

                    //----------------------------------------------------------
                    /************* lij is nonzero, x[i] is zero****************/
                    // x[i] = 0 then only perform IPGE update
                    // subtraction/division
                    //----------------------------------------------------------

                    if (sgn == 0)
                    {
                        // No previous pivot
                        if (jnew < 1)
                        {
                            // x[i] = 0 - lij*x[j]
                            SLIP_CHECK(SLIP_mpz_submul(x[i], L->x[m], x[j]));
                            h[i] = jnew;   // Entry is up to date
                        }

                        // Previous pivot exists
                        else
                        {
                            // x[i] = 0 - lij*x[j]
                            SLIP_CHECK(SLIP_mpz_submul(x[i], L->x[m], x[j]));

                            // x[i] = x[i] / rho[j-1]
                            SLIP_CHECK(SLIP_mpz_divexact(x[i], x[i],
                                rhos[jnew-1]));
                            h[i] = jnew;   // Entry is up to date
                        }
                    }

                    //----------------------------------------------------------
                    /************ Both lij and x[i] are nonzero****************/
                    // x[i] != 0 --> History & IPGE update on x[i]
                    //----------------------------------------------------------
                    else
                    {
                        // No previous pivot in this case
                        if (jnew < 1)
                        {
                            // x[i] = x[i]*rho[0]
                            SLIP_CHECK(SLIP_mpz_mul(x[i],x[i],rhos[0]));

                            // x[i] = x[i] - lij*xj
                            SLIP_CHECK(SLIP_mpz_submul(x[i], L->x[m], x[j]));
                            h[i] = jnew;   // Entry is now up to date
                        }
                        // There is a previous pivot
                        else
                        {
                            // History update if necessary
                            if (h[i] < jnew - 1)
                            {
                                // x[i] = x[i] * rho[j-1]
                                SLIP_CHECK(SLIP_mpz_mul(x[i], x[i],
                                    rhos[jnew-1]));
                                if (h[i] > -1)
                                {
                                    // x[i] = x[i] / rho[h[i]]
                                    SLIP_CHECK(SLIP_mpz_divexact(x[i], x[i],
                                        rhos[h[i]]));
                                }
                            }
                            // x[i] = x[i] * rho[j]
                            SLIP_CHECK(SLIP_mpz_mul(x[i],x[i],rhos[jnew]));
                            // x[i] = x[i] - lij*xj
                            SLIP_CHECK(SLIP_mpz_submul(x[i], L->x[m], x[j]));
                            // x[i] = x[i] / rho[j-1]
                            SLIP_CHECK(SLIP_mpz_divexact(x[i], x[i],
                                rhos[jnew-1]));
                            h[i] = jnew;   // Entry is up to date
                        }
                    }
                }
            }
        }
        else                               // Entries of L
        {
            //------------------------------------------------------------------
            // History update
            //------------------------------------------------------------------
            if (h[j] < k-1)
            {
                // x[j] = x[j] * rho[k-1]
                SLIP_CHECK(SLIP_mpz_mul(x[j], x[j], rhos[k-1]));
                if (h[j] > -1)
                {
                    // x[j] = x[j] / rho[h[j]]
                    SLIP_CHECK(SLIP_mpz_divexact(x[j], x[j], rhos[h[j]]));
                }
            }
        }
    }
    // Output the beginning of nonzero pattern
    *top_output = top;
    return SLIP_OK;
}

