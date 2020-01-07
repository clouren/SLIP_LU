#include "demos.h"

//------------------------------------------------------------------------------
// This example shows how to use multiple RHS vectors
//------------------------------------------------------------------------------

// usage:
// example4 > out
// out is file for output calculated result

int Ap[5] = {0, 3, 5, 8, 11};
int Ai[11]       = {0, 1, 2, 2, 3, 1, 2, 3, 0, 1,  2};
double Axnum[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};    // Numerator of x
double Axden[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};    // Denominator of x
double bxnum[4] = {17, 182, 61, 67};                      // Numerator of b
double bxden[4] = {15,  3,   6,  7};                      // Denominator of b
double bxnum2[4] = {8, 50, 25, 23};                       // Numerator of b2
double bxden2[4] = {15, 3, 6, 7};                         // Denominator of b2

#define FREE_WORKSPACE                          \
    SLIP_free(i);                               \
    SLIP_free(p);                               \
    SLIP_free(x);                               \
    SLIP_delete_double_mat(&b_doub, n, numRHS); \
    SLIP_delete_double_mat(&soln, n, numRHS);   \
    SLIP_delete_dense(&b);                      \
    SLIP_delete_LU_analysis(&S);                \
    SLIP_delete_sparse(&A);                     \
    SLIP_FREE(option);                          \
    SLIP_finalize( ) ;

int main (int argc, char **argv)
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
    int n = 4, nz = 11, numRHS = 2;
    double* x = NULL;
    double **b_doub = NULL;
    int* i = NULL;
    int* p = NULL;
    double** soln = NULL;
    SLIP_sparse* A = SLIP_create_sparse();
    SLIP_dense *b = SLIP_create_dense();
    SLIP_LU_analysis* S = SLIP_create_LU_analysis(n+1);
    SLIP_options* option = SLIP_create_default_options();
    if (!A || !b || !S || !option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }
    //--------------------------------------------------------------------------
    // Get matrix
    //--------------------------------------------------------------------------
    x = (double*) SLIP_malloc(nz* sizeof(double));
    b_doub = SLIP_create_double_mat(n, numRHS);
    i = (int*) SLIP_malloc(nz* sizeof(int));
    p = (int*) SLIP_malloc((n+1)* sizeof(int));
    soln = SLIP_create_double_mat(n, numRHS);;
    if (!x || !b_doub || !i || !p || !soln)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }
    for (int j = 0; j < n; j++)                           // Get b & p
    {
        p[j] = Ap[j];
        b_doub[j][0] = bxnum[j]/bxden[j];
        b_doub[j][1] = bxnum2[j]/bxden2[j];
    }
    p[n] = Ap[n];
    for (int j = 0; j < nz; j++)                          // Get i and x
    {
        i[j] = Ai[j];
        x[j] = Axnum[j]/Axden[j];
    }
    //--------------------------------------------------------------------------
    // Build A and b
    //--------------------------------------------------------------------------
    OK(SLIP_build_sparse_ccf_double(A, p, i, x, n, nz, option));
    OK(SLIP_build_dense_double(b, b_doub, n, numRHS, option));

    //--------------------------------------------------------------------------
    // Factorize
    //--------------------------------------------------------------------------

    clock_t start_sym = clock();

    // Symbolic analysis to obtain the column ordering
    OK(SLIP_LU_analyze(S, A, option));

    clock_t end_sym = clock();

    clock_t start_f = clock();

    // Solve the linear system using SLIP LU. The keyword double below indicates
    // that the solution will be given as a double**
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
