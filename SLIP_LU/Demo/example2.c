#include "demos.h"
// This example shows how to use SLIP LU within your code and read in a matrix
// stored in MM format. Also shows how to use SLIP with an mpq output

// usage:
// example2 mat_file rhs_file > out
// mat_file is the Matrix Market file containing the A matrix
// rhs_file is a list of entries for right hand side dense matrix
// if input file names are not specified, they are defaulted to
// ../ExampleMats/10teams.mat and ../ExampleMats/10teams.v, respectively.
// out is file for output calculated result


#define FREE_WORKSPACE                  \
    SLIP_delete_LU_analysis(&S);        \
    SLIP_delete_sparse(&A);             \
    SLIP_FREE(option);                  \
    SLIP_delete_dense(&b);              \
    SLIP_delete_mpq_mat(&x, n, numRHS); \
    SLIP_finalize();


int main (int argc, char **argv)
{
    //--------------------------------------------------------------------------
    // Prior to using SLIP LU, its environment must be initialized. This is done
    // by calling the SLIP_initialize() function. 
    //--------------------------------------------------------------------------
    SLIP_initialize();
    SLIP_info ok;
    int n = 0, numRHS = 0;

    //--------------------------------------------------------------------------
    // Get matrix and right hand side file names
    //--------------------------------------------------------------------------
    char *mat_name, *rhs_name;
    mat_name = "../ExampleMats/10teams_mat.txt";
    rhs_name = "../ExampleMats/10teams_v.txt";
    if (argc > 2)
    {
        mat_name = argv[1];
        rhs_name = argv[2];
    }

    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    mpq_t** x = NULL;
    SLIP_sparse *A = SLIP_create_sparse();
    SLIP_dense *b = SLIP_create_dense();
    SLIP_options* option = SLIP_create_default_options();
    SLIP_LU_analysis* S = NULL;
    if (!A || !b || !option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Allocate memory, read in A and b
    //--------------------------------------------------------------------------
    // Read in A
    FILE* mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    OK(SLIP_tripread(A, mat_file));
    fclose(mat_file);

    // Read in right hand side
    FILE* rhs_file = fopen(rhs_name,"r");
    if( rhs_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    OK(SLIP_read_dense(b, rhs_file));
    fclose(rhs_file);

    // Check if the size of A matches b
    if (A->m != b->m)
    {
        printf("%d %d \n", A->m,b->m);
        fprintf (stderr, "Error! Size of A and b do not match!\n");
        FREE_WORKSPACE;
        return 0;
    }
    n = A->m;                                             // Set n
    numRHS = b->n;
    S = SLIP_create_LU_analysis(n+1);
    x = SLIP_create_mpq_mat(n,numRHS);
    if (!S || !x)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Symbolic Ordering and Factorization
    //--------------------------------------------------------------------------

    clock_t start_sym = clock();

    // Symbolic analysis to obtain column permutation
    OK(SLIP_LU_analyze(S, A, option));

    clock_t end_sym = clock();

    clock_t start_f = clock();

    // Solve the linear system using SLIP LU. The keyword mpq below indicates
    // that the final solution vector x will be given as a mpq_t**
    OK(SLIP_solve_mpq(x, A, S, b, option));

    clock_t end_f = clock();

    double t_s = (double) (end_sym - start_sym) / CLOCKS_PER_SEC;
    double t_f = (double) (end_f - start_f) / CLOCKS_PER_SEC;

    printf ("\nSymbolic Analysis Time: %lf", t_s);
    printf ("\nSLIP LU Factor & Solve time: %lf\n", t_f);

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;

    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}
