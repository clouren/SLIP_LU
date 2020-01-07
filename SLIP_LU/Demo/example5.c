#include "demos.h"

// usage:
// example5 mat_file > out
// mat_file is the Matrix Market file containing the A matrix
// if input file name is not specified, it is defaulted to
// ../ExampleMats/10teams_mat.txt
// out is file for output calculated result


#define FREE_WORKSPACE                           \
    SLIP_delete_double_mat(&b_doub, n, numRHS);  \
    SLIP_delete_double_mat(&soln, n, numRHS);    \
    SLIP_delete_dense(&b);                       \
    SLIP_delete_LU_analysis(&S);                 \
    SLIP_delete_sparse(&A);                      \
    SLIP_FREE(option);                           \
    SLIP_finalize ( ) ;


int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // Prior to using SLIP LU, its environment must be initialized. This is done
    // by calling the SLIP_initialize() function. 
    //--------------------------------------------------------------------------
    SLIP_initialize();

    //--------------------------------------------------------------------------
    // Get matrix file names
    //--------------------------------------------------------------------------
    char *mat_name;
    if (argc > 1)
    {
        mat_name = argv[1];
    }
    else
    {
        mat_name = "../ExampleMats/10teams_mat.txt";
    }

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------
    SLIP_info ok;
    int n = 177, numRHS = 500;
    double **b_doub = NULL;
    double **soln = NULL;
    SLIP_sparse *A = SLIP_create_sparse();
    SLIP_dense *b = SLIP_create_dense();
    SLIP_LU_analysis *S = SLIP_create_LU_analysis(n+1);
    SLIP_options* option = SLIP_create_default_options();
    if (!A || !b || !S || !option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Generate random b_doub
    //--------------------------------------------------------------------------
    b_doub = SLIP_create_double_mat(n, numRHS);;
    soln = SLIP_create_double_mat(n, numRHS);
    if (!b_doub || !soln)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }
    int seed = 12;      // Just an arbitrary random seed
    srand (seed);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            b_doub[i][j] = rand();
        }
    }

    //--------------------------------------------------------------------------
    // Read in the matrix specified by mat_name and store it in A
    //--------------------------------------------------------------------------
    FILE* mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    OK(SLIP_tripread(A, mat_file));
    fclose(mat_file);

    //--------------------------------------------------------------------------
    // Build b
    //--------------------------------------------------------------------------
    OK(SLIP_build_dense_double(b, b_doub, n, numRHS, option));

    //--------------------------------------------------------------------------
    // Factorize
    //--------------------------------------------------------------------------

    clock_t start_sym = clock();

    // Symbolic analysis to obtain the column ordering of A
    OK(SLIP_LU_analyze(S, A, option));

    clock_t end_sym = clock();

    clock_t start_f = clock();

    // Solve the linear system using the SLIP LU factorization. The keyword double 
    // below indicates that the final solution will be returned as double**
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
