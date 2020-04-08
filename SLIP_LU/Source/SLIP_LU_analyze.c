//------------------------------------------------------------------------------
// SLIP_LU/SLIP_LU_analyze: symbolic ordering and analysis for sparse LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs the symbolic ordering for SLIP LU.
 * Currently, there are three options: user-defined order, COLAMD, or AMD.
 *
 * Input/output arguments:
 *
 * S:       Symbolic analysis struct. Undefined on input; contains column
 *          permutation and guesses on L and U nnz on output
 *
 * A:       Input matrix, unmodified on input/output
 *
 * option:  option->order tells the function which ordering scheme to use
 *
 */

// SLIP_LU_analyze creates the SLIP_LU_analysis object S.  Use
// SLIP_delete_LU_analysis to delete it.

#include "slip_internal.h"

SLIP_info SLIP_LU_analyze
(
    SLIP_LU_analysis** S_handle, // symbolic analysis (column perm. and nnz L,U)
    const SLIP_matrix *A,        // Input matrix
    const SLIP_options *option   // Control parameters, if NULL, use default
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // A can have any data type, but must be in sparse CSC format
    SLIP_REQUIRE_KIND (A, SLIP_CSC) ;

    if (!S_handle || A->n != A->m)
    {
        return SLIP_INCORRECT_INPUT;
    }
    (*S_handle) = NULL ;

    //--------------------------------------------------------------------------
    // allocate symbolic analysis object
    //--------------------------------------------------------------------------

    SLIP_LU_analysis *S = NULL ;
    int64_t i, n = A->n, anz = SLIP_matrix_nnz(A, option);
    // ALlocate memory for S
    S = (SLIP_LU_analysis*) SLIP_malloc(sizeof(SLIP_LU_analysis));
    if (S == NULL) {return SLIP_OUT_OF_MEMORY;}

    // Allocate memory for column permutation
    S->q = (int64_t*) SLIP_malloc((n+1) * sizeof(int64_t));
    if (S->q == NULL)
    {
        SLIP_FREE(S);
        return SLIP_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // No ordering is used. S->q is set to [0 ... n] and the number of nonzeros
    // in L and U is estimated to be 10 times the number of nonzeros in A. This
    // is a very crude estimate on the nnz(L) and nnz(U)
    //--------------------------------------------------------------------------

    SLIP_col_order order = SLIP_OPTION_ORDER (option) ;
    int pr = SLIP_OPTION_PRINT_LEVEL (option) ;

    if (order == SLIP_NO_ORDERING)
    {
        for (i = 0; i < n+1; i++)
        {
            S->q[i] = i;
        }
        // Guesses for number of L and U nonzeros
        S->lnz = S->unz = 10*anz;
    }

    //--------------------------------------------------------------------------
    // The AMD ordering is used. S->q is set to AMD's column ordering on
    // A+A'. The number of nonzeros in L and U is given as AMD's computed
    // number of nonzeros in the Cholesky factor L of A+A'
    //--------------------------------------------------------------------------
    else if (order == SLIP_AMD)
    {
        double Control [AMD_CONTROL];           // Declare AMD control
        amd_l_defaults (Control) ;              // Set AMD defaults
        double Info [AMD_INFO];
        amd_l_order(n, A->p, A->i, S->q, Control, Info); // Perform AMD
        S->lnz = S->unz = Info[AMD_LNZ];        // Guess for unz and lnz
        if (pr > 0)   // Output AMD info if desired
        {
            SLIP_PRINTF ("\n****Column Ordering Information****\n") ;
            amd_l_control (Control) ;
            amd_l_info (Info) ;
        }
    }

    //--------------------------------------------------------------------------
    // The COLAMD ordering is used. S->q is set as COLAMD's column ordering.
    // The number of nonzeros in L and U is set as 10 times the number of
    // nonzeros in A. This is a crude estimate.
    //--------------------------------------------------------------------------
    else
    {
        // Declared as per COLAMD documentation
        int64_t Alen = 2*anz + 6 *(n+1) + 6*(n+1) + n;
        int64_t* A2 = (int64_t*) SLIP_malloc(Alen* sizeof(int64_t));
        if (!A2)
        {
            // out of memory
            SLIP_LU_analysis_free (&S, option) ;
            return (SLIP_OUT_OF_MEMORY) ;
        }
        // Initialize S->q as per COLAMD documentation
        for (i = 0; i < n+1; i++)
        {
            S->q[i] = A->p[i];
        }
        // Initialize A2 per COLAMD documentation
        for (i = 0; i < anz; i++)
        {
            A2[i] = A->i[i];
        }
        int64_t stats [COLAMD_STATS];
        colamd_l (n, n, Alen, A2, S->q, (double *) NULL, stats);
        // Guess for lnz and unz
        S->lnz = S->unz = 10*anz;

        // Print stats if desired
        if (pr > 0)
        {
            SLIP_PRINTF ("\n****Column Ordering Information****\n") ;
            colamd_l_report (stats) ;
            SLIP_PRINTF ("\nEstimated L and U nonzeros: %" PRId64 "\n", S->lnz);
        }
        SLIP_FREE(A2);
    }

    //--------------------------------------------------------------------------
    // Make sure appropriate space is allocated. It's possible to return
    // guesses which exceed the dimension of L and U or guesses which are too
    // small for L U. In this case, this block of code ensures that the guesses
    // on nnz(L) and nnz(U) are at least n and no more than n*n.
    //--------------------------------------------------------------------------
    // Guess exceeds max number of nnz in A
    if (S->lnz > (double) n*n)
    {
        int64_t nnz = ceil(0.5*n*n);
        S->lnz = S->unz = nnz;
    }
    // If < n, first column of triangular solve may fail
    if (S->lnz < n)
    {
        S->lnz = S->lnz + n;
    }
    if (S->unz < n)
    {
        S->unz = S->unz + n;
    }

    //--------------------------------------------------------------------------
    // return result
    //--------------------------------------------------------------------------

    (*S_handle) = S ;
    return SLIP_OK;
}

