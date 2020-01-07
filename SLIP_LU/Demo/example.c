#include "demos.h"

/* This example shows how to use SLIP LU with a given input matrix and a double
   output */

// usage:
// example > out
// out is file for output calculated result

#define FREE_WORKSPACE                         \
    SLIP_free(i);                              \
    SLIP_free(j);                              \
    SLIP_free(x);                              \
    SLIP_delete_double_mat(&soln, n, 1);       \
    SLIP_delete_double_mat(&b_doub, n, 1);     \
    SLIP_delete_dense(&b);                     \
    SLIP_delete_LU_analysis(&S);               \
    SLIP_delete_sparse(&A);                    \
    SLIP_FREE(option);                         \
    SLIP_finalize() ;

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
    int n = 50, nz = 2500, num = 0;
    int *i = NULL;
    int *j = NULL;
    double *x = NULL;
    double **b_doub = NULL;
    double **soln = NULL;                       // Solution vector
    SLIP_sparse *A = SLIP_create_sparse();      // input matrix
    SLIP_dense *b = SLIP_create_dense();       // Right hand side vector
    SLIP_LU_analysis *S = SLIP_create_LU_analysis(n+1);// Column permutation
    SLIP_options *option = SLIP_create_default_options();
    if (!A || !b || !S || !option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }
    //--------------------------------------------------------------------------
    // Generate a random dense 50*50 matrix
    //--------------------------------------------------------------------------
    i = (int*) SLIP_malloc(nz* sizeof(int));
    j = (int*) SLIP_malloc(nz* sizeof(int));
    x = (double*) SLIP_malloc(nz* sizeof(double));
    b_doub = SLIP_create_double_mat(n,1);
    soln = SLIP_create_double_mat(n,1);
    if (!i || !j || !x || !b_doub || !soln)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    int seed = 10;
    srand(seed);
    for (int k = 0; k < n; k++)
    {
        b_doub[k][0] = rand();
        for (int p = 0; p < n; p++)
        {
            i[num] = k;
            j[num] = p;
            x[num] = rand();
            num+=1;
        }
    }

    //--------------------------------------------------------------------------
    // Build A and b
    //--------------------------------------------------------------------------
    OK(SLIP_build_sparse_trip_double(A, i, j, x, n, nz, option));
    OK(SLIP_build_dense_double(b, b_doub, n, 1, option));

    //--------------------------------------------------------------------------
    // Factorize & solve
    //--------------------------------------------------------------------------

    clock_t start_sym = clock();

    // Column permutation of A
    OK(SLIP_LU_analyze(S, A, option));

    clock_t end_sym = clock();

    clock_t start_f = clock();

    // Solve the linear system. The keyword double below indicates that the final 
    // solution vector will be output as a double**
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
