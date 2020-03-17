//------------------------------------------------------------------------------
// SLIP_LU/Demo/example.c: example main program for SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

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
    int32_t n = 50, nz = 2500, num = 0;
    int32_t *i = NULL;
    int32_t *j = NULL;
    double *x = NULL;
    double **b_doub = NULL;
    double **soln = NULL;                       // Solution vector
    SLIP_sparse *A = NULL ;                     // input matrix
    SLIP_dense *b = NULL ;                      // Right hand side vector
    SLIP_LU_analysis *S = NULL ;                // Column permutation
    SLIP_options *option = SLIP_create_default_options();
    if (!option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Generate a random dense 50*50 matrix
    //--------------------------------------------------------------------------

    i = (int32_t*) SLIP_malloc(nz* sizeof(int32_t));
    j = (int32_t*) SLIP_malloc(nz* sizeof(int32_t));
    x = (double*) SLIP_malloc(nz* sizeof(double));
    b_doub = SLIP_create_double_mat(n,1);
    soln = SLIP_create_double_mat(n,1);
    if (!i || !j || !x || !b_doub || !soln)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    int32_t seed = 10;
    srand(seed);
    for (int32_t k = 0; k < n; k++)
    {
        b_doub[k][0] = rand();
        for (int32_t p = 0; p < n; p++)
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

    OK (SLIP_build_sparse_trip_double (&A, i, j, x, n, nz, option));
    OK (SLIP_build_dense_double (&b, b_doub, n, 1, option));

    //--------------------------------------------------------------------------
    // analyze and solve
    //--------------------------------------------------------------------------

    clock_t start_sym = clock();

    // Column permutation of A
    OK(SLIP_LU_analyze(&S, A, option));

    clock_t end_sym = clock();

    clock_t start_f = clock();

    // Solve the linear system. The keyword double below indicates that the
    // final solution vector will be output as a double**.
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

