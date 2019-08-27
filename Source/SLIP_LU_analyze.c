# include "SLIP_LU_internal.h"

/*
 * Purpose: This function performs the symbolic ordering for SLIP LU. Currently,
 * there are four options: user defined order, COLAMD, AMD, or UMFPACK.
 */
#define SLIP_FREE_WORKSPACE   \
    SLIP_FREE(p);             \
    SLIP_FREE(q);             \
    SLIP_FREE(Ax);


SLIP_info SLIP_LU_analyze
(
    SLIP_LU_analysis* S,  // symbolic analysis (column permutation and nnz L,U)
    SLIP_sparse* A,       // Input matrix
    SLIP_dense *b,        // right hand side vectors (UMFPACK only)
    SLIP_options* option  // Control parameters
)
{
    // do not check for the availability of b here, which will be checked in 
    // UMFPACK. For other ordering methods, b can be NULL
    if (!S || !A || !(A->i) || !(A->x) || !(A->p) || !option || A->n != A->m)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t n = A->n, nz = A->nz, i, k;
    SLIP_info ok;
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
        S->lnz = 10*nz; S->unz = S->lnz;
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
        S->lnz = Info[AMD_LNZ]; S->unz = S->lnz;// Guess for unz and lnz
        if (option->print_level > 0)            // Output AMD info if desired
        {
            printf("\n****Column Ordering Information****\n");
            amd_control(Control);
            amd_info(Info);
        }
    }

    //--------------------------------------------------------------------------
    // UMFPACK
    //--------------------------------------------------------------------------
    else if (option->order == SLIP_UMFPACK)
    {
        int32_t do_recip;
        // Must override pivoting scheme when using UMFPACK
        option->pivot = SLIP_DIAGONAL;
        int32_t* p = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
        int32_t* q = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
        // Numeric values of Ax in double precision
        double* Ax = (double*) SLIP_malloc(nz* sizeof(double));
        if (!p || !q || !Ax)
        {
	    SLIP_FREE_WORKSPACE;
            return SLIP_OUT_OF_MEMORY;
        }
        for (i = 0; i < nz; i++)
        {
             SLIP_CHECK(slip_mpz_get_d(&(Ax[i]), A->x[i]));
        }
        double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
        void* Symbolic, *Numeric;

        //----------------------------------------------------------------------
        // Do UMFPACK
        //----------------------------------------------------------------------
        // Set defaults
        umfpack_di_defaults(Control);
        // Unsymmetric ordering strategy
        Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
        // Pivot tolerance of 10^-6
        Control[UMFPACK_PIVOT_TOLERANCE] = 0.000001;
        // Evaluates AMD, COLAMD, Metis and Nested Dissection
        Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
        // UMFPACK Symbolic
        umfpack_di_symbolic(n, n, A->p, A->i, Ax, &Symbolic, Control, Info);
        // UMFPACK Numeric
        umfpack_di_numeric(A->p, A->i, Ax, Symbolic, &Numeric, Control, Info);
        // Get P and Q of UMFPACK
        umfpack_di_get_numeric(NULL, NULL, NULL, NULL, NULL, NULL, p,
            q, NULL, &do_recip, NULL, Numeric);

        // Free memory of symbolic and numeric
        umfpack_di_free_symbolic (&Symbolic);
        umfpack_di_free_numeric (&Numeric);

        // Set S->q as UMFPACK q
        for (k = 0; k < n; k++)
        {
            S->q[k] = q[k];
        }
        // UMFPACK lnz < SLIP lnz usually
        S->lnz = 2 * Info[UMFPACK_LNZ];
        // UMFPACK unz < SLIP unz usually
        S->unz = 2 * Info[UMFPACK_UNZ];
        SLIP_CHECK(slip_UMFPACK_permute(b, A, p));
	SLIP_FREE_WORKSPACE;
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
        S->lnz = 10*A->nz; S->unz = S->lnz;

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
        S->lnz = nnz; S->unz = nnz;
    }
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
#undef SLIP_FREE_WORKSPACE
