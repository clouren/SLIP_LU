//------------------------------------------------------------------------------
// SLIP_LU/Demo/example3.c: example main program for SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "demos.h"

//------------------------------------------------------------------------------
// This example uses an CSC input in mpfr precision/
//------------------------------------------------------------------------------

// usage:
// example3 > out
// out is file for output calculated result

int32_t Ap[5] = {0, 3, 5, 8, 11};
int32_t Ai[11]       = {0, 1, 2, 2, 3, 1, 2, 3, 0, 1,  2};
double Axnum[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};    // Numerator of x
double Axden[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};    // Denominator of x
double bxnum[4] = {17, 182, 61, 67};                      // Numerator of b
double bxden[4] = {15,  3,   6,  7};                      // Denominator of b

#define FREE_WORKSPACE                       \
    SLIP_free(i);                            \
    SLIP_free(p);                            \
    SLIP_delete_double_mat(&soln, n, 1);     \
    SLIP_delete_mpfr_mat(&b_mpfr, n, 1);     \
    SLIP_delete_mpfr_array(&x_mpfr, nz);     \
    SLIP_delete_dense(&b);                   \
    SLIP_delete_LU_analysis(&S);             \
    SLIP_delete_sparse(&A);                  \
    SLIP_FREE(option);                       \
    SLIP_finalize();

int main (void)
{

    //--------------------------------------------------------------------------
    // Prior to using SLIP LU, its environment must be initialized. This is done
    // by calling the SLIP_initialize() function.
    //--------------------------------------------------------------------------

    SLIP_initialize();

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------

    SLIP_info ok;
    int32_t n = 4, nz = 11, j;
    mpfr_t ** b_mpfr = NULL;
    int32_t* i = NULL;
    int32_t* p = NULL;
    mpfr_t* x_mpfr = NULL;
    double** soln = NULL;
    SLIP_sparse* A = NULL ;
    SLIP_dense *b = NULL ;
    SLIP_LU_analysis *S = NULL ;
    SLIP_options* option = SLIP_create_default_options();
    if (!option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Get matrix
    //--------------------------------------------------------------------------

    b_mpfr = SLIP_create_mpfr_mat(n, 1, option);
    i = (int32_t*) SLIP_malloc(nz* sizeof(int32_t));
    p = (int32_t*) SLIP_malloc((n+1)* sizeof(int32_t));
    x_mpfr = SLIP_create_mpfr_array(nz, option);
    soln = SLIP_create_double_mat(n,1);
    if (!b_mpfr || !i || !p || !x_mpfr || !soln)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    for (j = 0; j < n; j++)                               // Get p & b_mpfr
    {
        p[j] = Ap[j];
        mpfr_set_d(b_mpfr[j][0], bxnum[j], MPFR_RNDN);
        mpfr_div_d(b_mpfr[j][0], b_mpfr[j][0], bxden[j], MPFR_RNDN);
    }
    p[n] = Ap[n];                                        // Finalize p
    for (j = 0; j < nz; j++)                             // Get i and x_mpfr
    {
        i[j] = Ai[j];
        mpfr_set_d(x_mpfr[j], Axnum[j], MPFR_RNDN);
        mpfr_div_d(x_mpfr[j], x_mpfr[j], Axden[j], MPFR_RNDN);
    }

    //--------------------------------------------------------------------------
    // Build A and b
    //--------------------------------------------------------------------------

    OK(SLIP_build_sparse_csc_mpfr(&A, p, i, x_mpfr, n, nz, option));
    OK(SLIP_build_dense_mpfr(&b, b_mpfr, n, 1, option));

    //--------------------------------------------------------------------------
    // analyze
    //--------------------------------------------------------------------------

    clock_t start_sym = clock();

    // Symbolic analysis
    OK(SLIP_LU_analyze(&S, A, option));

    clock_t end_sym = clock();

    clock_t start_f = clock();

    //--------------------------------------------------------------------------
    // solve
    //--------------------------------------------------------------------------

    // Solve the linear system using SLIP LU. The keyword double below indicates
    // that the final solution vector will be given as double**
    OK(SLIP_solve_double(soln, A, S, b, option));

    clock_t end_f = clock();

    double t_s = (double) (end_sym - start_sym) / CLOCKS_PER_SEC;
    double t_f = (double) (end_f - start_f) / CLOCKS_PER_SEC;

    printf("\nSymbolic Analysis Time: %lf", t_s);
    printf("\nSLIP LU Factor & Solve time: %lf\n", t_f);

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;
    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}

