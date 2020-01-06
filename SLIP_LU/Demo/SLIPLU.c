#include "demos.h"

/* This program will exactly solve the sparse linear system Ax = b by performing
 * the SLIP LU factorization. Please refer to README.txt for information on how
 * to properly use this code
 */

// usage:
// SLIPLU Followed by the listed args:
//
// help. e.g., SLIPLU help, which indicates SLIPLU to print to guideline
// for using this function.
//
// f (or file) Filename. e.g., SLIPLU f MATRIX_NAME RHS_NAME, which indicates
// SLIPLU will read matrix from MATRIX_NAME and right hand side from RHS_NAME.
// The matrix must be stored in Matrix Market format. Please refer to
// http://math.nist.gov/MatrixMarket/formats.html for information on
// Matrix Market format.
// The right hand side vector must be stored as a dense vector.
//
// c (or check). e.g., SLIPLU c, which indicates SLIPLU will check
// the correctness of the solution via A*x == b.
// WARNING: This is slow and should only be used for verification
//
// p (or piv) Pivot_param. e.g., SLIPLU p 0, which inidcates SLIPLU will use
// smallest pivot for pivot scheme. Other available options are listed
// as follows:
//        0: Smallest pivot: Default and recommended
//        1: Diagonal pivoting
//        2: First nonzero per column chosen as pivot
//        3: Diagonal pivoting with tolerance for smallest pivot
//        4: Diagonal pivoting with tolerance for largest pivot
//        5: Largest pivot
//
// q (or col) Column_order_param. e.g., SLIPLU q 0, which indicates SLIPLU
// will use COLAMD for column ordering. Other available options are:
//        0: None: Not recommended for sparse matrices
//        1: COLAMD: Default
//        2: AMD
//
// t (or tol) tolerance_param. e.g., SLIPLU t 1e-10, which indicates SLIPLU
// will use 1e-10 as the tolerance for pivot scheme 3 and 4 mentioned above.
// Therefore, it is only necessary if pivot scheme 3 or 4 is used.
//
// o2 (or out2). e.g., SLIPLU o2 1, which indicates SLIPLU will output the
// errors and warnings during the process. Other available options are:
//        0: print nothing
//        1: just errors and warnings: Default
//        2: terse, with basic stats from COLAMD/AMD and SLIP and solution
//
// o (or out) output_param. e.g., SLIPLU o 1, which indicates SLIPLU will
// output the solution as full precision rational arithmetic stdout or file
// specified by the stdout redirection. It is noted that solution will be
// output only when o2 >= 2. Other available options are:
//        1: The solution will be output in full precision rational arithmetic
//        2: The solution will be output in double precision
//        3: The solution will be output in some user specified precision.
//           This must be input as o 3 USER_PRECISION. Precision must be
//           an integer, e.g., SLIPLU o 3 128.
//
// If none of the above args is given, they are set to the following default:
//
//  mat_name = "../ExampleMats/10teams_mat.txt"
//  rhs_name = "../ExampleMats/10teams_v.txt"
//  no solution checking
//  p = 0, i.e., using smallest pivot
//  q = 1, i.e., using COLAMD
//  t = 0.1, not being using since p != 3 or 4
//  o2 = 0, i.e., nothing will be printed
//  o = 1, not being using since o2 < 2


#define FREE_WORKSPACE                           \
    SLIP_delete_sparse(&A);                      \
    SLIP_delete_sparse(&L);                      \
    SLIP_delete_sparse(&U);                      \
    SLIP_delete_dense(&b);                       \
    SLIP_delete_mpz_array(&rhos,nrows);          \
    SLIP_delete_mpq_mat(&x,nrows,numRHS);        \
    SLIP_delete_double_mat(&x_doub,nrows,numRHS);\
    SLIP_delete_mpfr_mat(&x_mpfr,nrows,numRHS);  \
    SLIP_free(pinv);                             \
    SLIP_delete_LU_analysis(&S);                 \
    SLIP_FREE(option);                           \
    SLIP_finalize( ) ;

int main( int argc, char* argv[])
{
    SLIP_initialize();
    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int nrows = 0, numRHS = 0;
    int rat = 1;
    SLIP_info ok, check = SLIP_OK ;
    char *mat_name, *rhs_name;
    mat_name = "../ExampleMats/10teams_mat.txt";// Set demo matrix and RHS name
    rhs_name = "../ExampleMats/10teams_v.txt";

    mpz_t* rhos = NULL;
    mpq_t** x = NULL;
    int* pinv = NULL;
    SLIP_LU_analysis* S = NULL;
    double **x_doub = NULL;
    mpfr_t **x_mpfr = NULL;
    SLIP_sparse *A = SLIP_create_sparse();
    SLIP_sparse *L = SLIP_create_sparse();
    SLIP_sparse *U = SLIP_create_sparse();
    SLIP_dense *b  = SLIP_create_dense();
    SLIP_options *option = SLIP_create_default_options();
    if (!A || !L || !U || !b || !option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        SLIP_finalize();
        return 0;
    }
    // Process the command line
    OK(SLIP_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));

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
        fprintf (stderr, "Error! Size of A and b do not match!\n");
        FREE_WORKSPACE;
        return 0;
    }
    nrows = A->m;
    numRHS = b->n;

    rhos = SLIP_create_mpz_array(nrows);
    pinv = (int*) SLIP_malloc(nrows* sizeof(int));
    x = SLIP_create_mpq_mat(nrows, numRHS);
    S = SLIP_create_LU_analysis(nrows+1);
    if (!rhos || !pinv || !x || !S)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Perform Column ordering
    //--------------------------------------------------------------------------

    clock_t start_col = clock();

    // Column ordering using either AMD, COLAMD or nothing
    OK(SLIP_LU_analyze(S, A, option));
    if (option->print_level > 0)
    {
        SLIP_print_options(option);
    }

    clock_t end_col = clock();
    
    //--------------------------------------------------------------------------
    // SLIP LU Factorization
    //--------------------------------------------------------------------------
    clock_t start_factor = clock();

    OK(SLIP_LU_factorize(L, U, A, S, rhos, pinv, option));

    clock_t end_factor = clock();

    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    clock_t start_solve = clock();

    // Solve LDU x = b
    OK(SLIP_LU_solve(x, b, rhos, L, U, pinv));

    clock_t end_solve = clock();

    //--------------------------------------------------------------------------
    // Soln verification
    //--------------------------------------------------------------------------
    // x = Q x
    OK(SLIP_permute_x(x, nrows, numRHS, S));
    OK(SLIP_check_solution(A, x, b));

    // TODO
    check = ok;       // ok is assigned as the status of SLIP_check_solution
    if (check == SLIP_OK)
    {
        printf ("Solution is verified to be exact.\n") ;
    }
    else
    {
        printf ("ERROR! Solution is wrong.\n") ;
    }
    OK(SLIP_scale_x(x, A, b));

    //--------------------------------------------------------------------------
    // Output
    //--------------------------------------------------------------------------
    if (rat == 1)
    {
        OK(SLIP_print_stats_mpq(stdout, x, nrows, numRHS, check, option));
    }
    else if (rat == 2)
    {
        x_doub = SLIP_create_double_mat(nrows, numRHS);
        if (x_doub == NULL) {OK(SLIP_OUT_OF_MEMORY);}
        OK(SLIP_get_double_soln(x_doub, x, nrows, numRHS));
        OK(SLIP_print_stats_double(stdout, x_doub, nrows, numRHS, check,
            option));
    }
    else
    {
        x_mpfr = SLIP_create_mpfr_mat(nrows, numRHS, option);
        if (x_mpfr == NULL) {OK(SLIP_OUT_OF_MEMORY);}
        OK(SLIP_get_mpfr_soln(x_mpfr, x, nrows, numRHS, option));
        OK(SLIP_print_stats_mpfr(stdout, x_mpfr, nrows, numRHS, check, option));
    }

    double t_sym = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_factor = (double) (end_factor - start_factor) / CLOCKS_PER_SEC;
    double t_solve =  (double) (end_solve - start_solve) / CLOCKS_PER_SEC;

    printf("\nNumber of L+U nonzeros: \t\t%d",
        (L->nz) + (U->nz) - (L->m));
    printf("\nSymbolic analysis time: \t\t%lf", t_sym);
    printf("\nSLIP LU Factorization time: \t\t%lf", t_factor);
    printf("\nFB Substitution time: \t\t\t%lf\n\n", t_solve);

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}
