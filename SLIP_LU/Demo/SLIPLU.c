#include "demos.h"

/* This program will exactly solve the sparse linear system Ax = b by performing
 * the SLIP LU factorization. This is intended to be a demonstration of the 
 * "advanced interface" of SLIP LU. Please refer to README.txt for information 
 * on how to properly use this code
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
// q (or col) Column_order_param. e.g., SLIPLU q 1, which indicates SLIPLU
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
    
    //--------------------------------------------------------------------------
    // Prior to using SLIP LU, its environment must be initialized. This is done
    // by calling the SLIP_initialize() function. 
    //--------------------------------------------------------------------------
    
    SLIP_initialize();
    
    //--------------------------------------------------------------------------
    // We first initialize the default parameters. These parameters are modified
    // either via command line arguments or when reading in data. The important 
    // initializations are in the second block below, where we initialize the following:
    //
    //  rhos: Sequence of pivots used in LU
    //
    //  x: The final solution vector of the system (notice it is rational mpq_t)
    //  pinv: inverse row permutation used for LU factorization
    //
    //  S: Analysis struct that contains the column ordering used
    //
    //  x_doub: Final solution to the linear system stored as double precision (note 
    //          that we declare it here to enumerate all possibilities).
    //
    //  x_mpfr: Final solution to the linear system stored as double precision (note 
    //          that we declare it here to enumerate all possibilities).
    //
    //  A, L, U: Sparse integer matrices. This is the default struct used within
    //           SLIP LU to perform most internal routines. Sparse matrices are
    //           always initialized with a call to SLIP_create_sparse(); Note that
    //           the numeric entries within these matrices are integral mpz_t data
    //           types. If the input matrix contains non integral entries, it must
    //           be appropriately scaled using one of the SLIP_build_sparse_* functions.
    //
    //  b: Dense right hand side vector(s). Currently SLIP LU assumes that the RHS
    //     is always dense. All dense matrices must be initialized with the function
    //     SLIP_create_dense(); Note that dense matrices are subject to the same caveat
    //     as the sparse matrices, namely that their entries are assumed to be integral
    //     thus if they are not they must be appropraite scaled using the appropriate 
    //     SLIP_build_dense_* function.
    //
    //  option: The SLIP_options struct contains various command parameters for 
    //          the factorization which are highlighed both above and in the user 
    //          guide. Much like a sparse and dense matrix, it must be created 
    //          with a call to SLIP_create_default_options();
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
    
    //--------------------------------------------------------------------------
    // After initializing memory, we process the command line for this function.
    // Such a step is optional, a user can also manually set these parameters. 
    // For example, if one wished to use the AMD ordering, they can just set
    // option->order = SLIP_AMD.
    //--------------------------------------------------------------------------
    
    OK(SLIP_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));

    //--------------------------------------------------------------------------
    // In this demo file, we now read in the A and b matrices from external files.
    // Please refer to the example*.c files or the user for other methods of 
    // creating the input matrix
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
    if (A->n != b->m)
    {
        fprintf (stderr, "Error! Size of A and b do not match!\n");
        FREE_WORKSPACE;
        return 0;
    }
    nrows = A->m;
    numRHS = b->n;

    //--------------------------------------------------------------------------
    // Now that we have read in the input matrix, we allocate memory for the pivots,
    // inverse row permutation, solution to the system, and the analysis struct.
    //--------------------------------------------------------------------------
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
    // We now perform symbolic analysis by getting the column preordering of
    // the matrix A. This is done via the SLIP_LU_analyze function. The output
    // of this function is a column permutation Q where we factor the matrix AQ.
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
    // Now we perform the SLIP LU factorization to obtain matrices L and U and a
    // row permutation P such that PAQ = LDU. Note that the D matrix is never
    // explicitly constructed or used.
    //
    // Note that in the simple interface of SLIP LU (shown in the example*.c files)
    // it is not necessary to utilize the L and U matrices nor split the solution 
    // process into factorization and solve.
    //--------------------------------------------------------------------------
    clock_t start_factor = clock();

    OK(SLIP_LU_factorize(L, U, A, S, rhos, pinv, option));

    clock_t end_factor = clock();

    //--------------------------------------------------------------------------
    // We now solve the system Ax=b using the L and U factors computed above. 
    // After we obtain the solution, we permute it with respect to the column 
    // permutation (x_final = Q x).
    //--------------------------------------------------------------------------
    clock_t start_solve = clock();

    // Solve LDU x = b
    OK(SLIP_LU_solve(x, b, rhos, L, U, pinv));

    clock_t end_solve = clock();

    OK(SLIP_permute_x(x, nrows, numRHS, S));
    
    //--------------------------------------------------------------------------
    // SLIP LU has an optional check step which can verify that the solution vector
    // x satisfies Ax=b in perfect precision. 
    //
    // Note that this is entirely OPTIONAL and NOT NECESSARY. The solution returned
    // is guaranteed to be exact. Also, note that this function can be quite time
    // consuming; thus it is not recommended to be used in general.
    // 
    // This function returns the status SLIP_OK if it is successfully verified to
    // be correct. 
    //--------------------------------------------------------------------------
    
    check = SLIP_check_solution(A, x, b);

    if (check == SLIP_OK)
    {
        printf ("Solution is verified to be exact.\n") ;
    }
    else if (check == SLIP_INCORRECT)
    {
        // This should never happen.
        printf ("ERROR! Solution is wrong.\n") ;
    }
    else
    {
        // Out of memory or bad input.
        FREE_WORKSPACE;
        return 0;
    }
            
    //--------------------------------------------------------------------------
    // The x vector is now scaled. This consists of accounting for any scaling
    // done to A and b to make these entries integral. 
    //--------------------------------------------------------------------------
    OK(SLIP_scale_x(x, A, b));

    //--------------------------------------------------------------------------
    // Output timing statistics. This also prints to a file if desired.
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
