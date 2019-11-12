//------------------------------------------------------------------------------
// SLIP_LU/SLIP_LU_analyze: symbolic ordering and analysis for sparse LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/*
 * Purpose: This function performs the symbolic ordering for SLIP LU. Currently,
 * there are three options: user-defined order, COLAMD, or AMD.
 */

#define SLIP_FREE_WORKSPACE   \
    SLIP_FREE(p);             \
    SLIP_FREE(q);             \
    SLIP_FREE(Ax);

# include "SLIP_LU_internal.h"

SLIP_info SLIP_LU_analyze
(
    SLIP_LU_analysis* S,  // symbolic analysis (column permutation and nnz L,U)
    SLIP_sparse* A,       // Input matrix
    SLIP_options* option  // Control parameters
)
{
    if (!S || !A || !(A->i) || !(A->x) || !(A->p) || !option || A->n != A->m)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t n = A->n, nz = A->nz, i;
    // Print info if needed
    if (option->print_level > 0)
    {
        slip_lu_info();
    }

    //--------------------------------------------------------------------------
    // No ordering
    //--------------------------------------------------------------------------
    if (option->order == SLIP_NO_ORDERING)
    {
        for (i = 0; i < n+1; i++)
        {
            S->q[i] = i;
        }
        // Guesses for number of L and U nonzeros
        S->lnz = S->unz = 10*nz;
    }

    //--------------------------------------------------------------------------
    // AMD
    //--------------------------------------------------------------------------
    else if (option->order == SLIP_AMD)
    {
        double Control [AMD_CONTROL];           // Declare AMD control
        amd_defaults(Control);                  // Set AMD defaults
        double Info [AMD_INFO];
        amd_order(n, A->p, A->i, S->q, Control, Info); // Perform AMD
        S->lnz = S->unz = Info[AMD_LNZ];        // Guess for unz and lnz
        if (option->print_level > 0)            // Output AMD info if desired
        {
            printf("\n****Column Ordering Information****\n");
            amd_control(Control);
            amd_info(Info);
        }
    }

    //--------------------------------------------------------------------------
    // COLAMD
    //--------------------------------------------------------------------------
    else
    {
        // Declared as per COLAMD documentation
        int32_t Alen = 2*A->nz + 6 *(n+1) + 6*(n+1) + n;
        int32_t* A2 = (int32_t*) SLIP_malloc(Alen* sizeof(int32_t));
        if (!A2) {return SLIP_OUT_OF_MEMORY;}
        // Initialize S->q as per COLAMD documentation
        for (i = 0; i < n+1; i++)
        {
            S->q[i] = A->p[i];
        }
        // Initialize A2 per COLAMD documentation
        for (i = 0; i < nz; i++)
        {
            A2[i] = A->i[i];
        }
        int32_t stats [COLAMD_STATS];
        colamd(n, n, Alen, A2, S->q, (double *) NULL, stats);
        // Guess for lnz and unz
        S->lnz = S->unz = 10*A->nz;

        // Print stats if desired
        if (option->print_level > 0)
        {
            printf("\n****Column Ordering Information****\n");
            colamd_report(stats);
            printf("\nEstimated L and U nonzeros: %d", S->lnz);
        }
        SLIP_FREE(A2);
    }

    //--------------------------------------------------------------------------
    // Make sure too much space isn't allocated
    //--------------------------------------------------------------------------
    // Guess exceeds max number of nnz in A
    if (S->lnz > (double) n*n)
    {
        int32_t nnz = ceil(0.5*n*n);
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

    return SLIP_OK;
}

