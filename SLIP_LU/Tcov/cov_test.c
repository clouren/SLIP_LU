/*
 * When the test is run without input argument, brutal test is used and simple
 * test otherwise. Read the following for detailed instruction and information
 *
 *
 * For simple test, the test needs to be run with command
 * ./cov_test Ab_type rat N list1[0] ... list1[N-1] M list2[0] ... list2[M-1]
 * Ab_type: type of Matrix A and vector b: 0 mpz, 1 double, 2 int, 3 mpq, 4 mpfr
 *                                         5 for miscellaneous test
 * rat: type of solution: 1: full precision rational arithmetic,
 *                        2: double, 3:  user specified precision.
 * N and list1 specify the test list for slip_gmp_ntrials (in SLIP_gmp.h)
 * M and list2 specify the test list for malloc_count (in tcov_malloc_test.h)
 * N, list1, M, list2 are optional, but N and list1 are required when M and
 * list2 is wanted
 *
 * For brutal test, the test is run with command
 * ./cov_test
 * the test will run through all cases
 * (specifically, [Ab_type rat]={[0 1], [1 1], [2 2], [3 3], [4 3], [5 3]})
 * each case run from malloc_count = 0 to a number that can guarantee
 * malloc_count > 0 when the case finishes
 */

/* For simple test ONLY!!
 * uncomment to show the input lists for slip_gmp_ntrials and malloc_count
 */
// #define SLIP_TCOV_SHOW_LIST

/* This program will exactly solve the sparse linear system Ax = b by performing
 * the SLIP LU factorization. Please refer to README.txt for information on how
 * to properly use this code
 */

#define SLIP_FREE_WORKSPACE                      \
{                                                \
    SLIP_delete_LU_analysis(&S);                 \
    SLIP_delete_sparse(&A);                      \
    SLIP_FREE(option);                           \
    SLIP_delete_dense(&b);                       \
    SLIP_delete_mpz_mat(&B_mpz, n, numRHS);      \
    SLIP_delete_mpz_array(&Ax_mpz,nz);           \
    SLIP_FREE(Ax_doub);                          \
    SLIP_delete_double_mat(&B_doub, n, numRHS);  \
    SLIP_delete_int_mat(&B_int, n, numRHS);      \
    SLIP_FREE(Ax_int);                           \
    SLIP_delete_mpq_mat(&B_mpq , n, numRHS);     \
    SLIP_delete_mpq_array(&Ax_mpq ,nz);          \
    SLIP_delete_mpfr_mat(&B_mpfr, n, numRHS);    \
    SLIP_delete_mpfr_array(&Ax_mpfr,nz);         \
    SLIP_delete_mpq_mat(&sol_mpq, n, numRHS);    \
    SLIP_delete_double_mat(&sol_doub, n, numRHS);\
    SLIP_delete_mpfr_mat(&sol_mpfr, n, numRHS);  \
    SLIP_free(x_doub);                           \
    SLIP_delete_mpz_array(&x_mpz, nz);           \
    SLIP_delete_mpq_array(&x_mpq, nz);           \
    SLIP_delete_mpfr_array(&x_mpfr, nz);         \
    SLIP_delete_sparse(&M);                      \
    if (mat_file != NULL) {fclose(mat_file);}    \
    SLIP_finalize() ;                            \
}

#include "demos.h"

#define TEST_CHECK(method)                       \
{                                                \
    ok = method;                                 \
    if (ok!=SLIP_OK)                             \
    {                                            \
        SLIP_PRINT_OK(ok);                       \
        SLIP_FREE_WORKSPACE;                     \
        continue;                                \
    }                                            \
}

#define TEST_CHECK_FAILURE(method)               \
{                                                \
    ok = method;                                 \
    if (ok!=SLIP_INCORRECT_INPUT && ok!=SLIP_SINGULAR)\
    {                                            \
        SLIP_PRINT_OK(ok);                       \
        SLIP_FREE_WORKSPACE;                     \
        continue;                                \
    }                                            \
    else                                         \
    {                                            \
        printf("Expected failure at line %d\n", __LINE__);\
    }                                            \
}

#define CLEAR_SLIP_MAT_A                         \
{                                                \
    if (A != NULL)                               \
    {                                            \
        SLIP_delete_mpz_array(&(A->x), A->nzmax);\
        SLIP_FREE (A->i);                        \
        SLIP_FREE (A->p);                        \
        /*SLIP_mpq_set_ui(A->scale, 1, 1);*/     \
        A -> n     = 0;                          \
        A -> m     = 0;                          \
        A -> nzmax = 0;                          \
        A -> nz    = 0;                          \
    }                                            \
}

#define CLEAR_SLIP_MAT_M                         \
{                                                \
    if (M != NULL)                               \
    {                                            \
        SLIP_delete_mpz_array(&(M->x), M->nzmax);\
        SLIP_FREE (M->i);                        \
        SLIP_FREE (M->p);                        \
        /*SLIP_mpq_set_ui(M->scale, 1, 1);*/     \
        M -> n     = 0;                          \
        M -> m     = 0;                          \
        M -> nzmax = 0;                          \
        M -> nz    = 0;                          \
    }                                            \
}

#define CLEAR_SLIP_MAT_B                         \
{                                                \
    if (b != NULL)                               \
    {                                            \
        SLIP_delete_mpz_mat(&(b->x), b->m, b->n);\
        /*SLIP_mpq_set_ui(b->scale, 1, 1);*/     \
        b -> n     = 0;                          \
        b -> m     = 0;                          \
    }                                            \
}

int Ap[5] = {0, 3, 5, 8, 11};
int Ai[11]       = {0, 1, 2, 2, 3, 1, 2, 3, 0, 1,  2};
double Axnum[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};  // Numerator of x
double Axden[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};  // Denominator of x
double bxnum[4] = {170, 1820, 61, 670};                // Numerator of b
double bxden[4] = {15,  3,   6,  7};                    // Denominator of b
int Axnum3[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};    // Numerator of x
int Axden3[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};    // Denominator of x
int bxnum3[4] = {17, 182, 61, 67};                      // Numerator of b
int bxden3[4] = {15,  3,   6,  7};                      // Denominator of b

int main( int argc, char* argv[])
{
    bool IS_SIMPLE_TEST = true;
    int rat_list[6] = {1, 1, 2, 3, 3, 3};  // only used in brutal test
    int Ab_type = 0, rat = 1, NUM_OF_TRIALS = 0, NUM_OF_MALLOC_T = 0;
    int64_t *gmp_ntrial_list=NULL;         // only used in simple test
    int *malloc_trials_list=NULL;          // only used in simple test
    if (argc == 1)                         // brutal test
    {
        IS_SIMPLE_TEST = false;
        NUM_OF_TRIALS = 6;
    }
    else if (argc == 2)                     // incorrect input
    {
        fprintf(stderr, "incorrect input\n");
        return 0;
    }
    else                                   // simple test
    {
        IS_SIMPLE_TEST = true;

        int arg_count = 0;
        //type of Matrix A and vector b: 0 mpz, 1 double, 2 int, 3 mpq, 4 mpfr
        Ab_type = atoi(argv[++arg_count]);
        //type of solution: 1: full precision rational arithmetic,
        //                  2: double, 3:  user specified precision.
        rat=atoi(argv[++arg_count]);
        if (!argv[++arg_count])
        {
            NUM_OF_TRIALS=1;
            gmp_ntrial_list=SLIP_malloc(NUM_OF_TRIALS* sizeof(int64_t));
            gmp_ntrial_list[0]=-1;
            arg_count--;
        }
        else
        {
            NUM_OF_TRIALS=atoi(argv[arg_count]);
            gmp_ntrial_list=SLIP_malloc(NUM_OF_TRIALS* sizeof(int64_t));
            for (int k=0; k<NUM_OF_TRIALS; k++)
            {
                if (argv[++arg_count])
                {
                    gmp_ntrial_list[k]=atoi(argv[arg_count]);
                }
                else
                {
                    fprintf(stderr, "WARNING: MISSING gmp trial\n");
                    NUM_OF_TRIALS=1;
                    gmp_ntrial_list[0]=-1;
                    arg_count--;
                }
            }
        }
        if (!argv[++arg_count])
        {
            NUM_OF_MALLOC_T=1;
            malloc_trials_list=SLIP_malloc(NUM_OF_MALLOC_T* sizeof(int));
            malloc_trials_list[0]=1000;//INT_MAX;
        }
        else
        {
            NUM_OF_MALLOC_T=atoi(argv[arg_count]);
            malloc_trials_list=SLIP_malloc(NUM_OF_MALLOC_T* sizeof(int));
            for (int k=0; k<NUM_OF_MALLOC_T; k++)
            {
                if (argv[++arg_count])
                {
                    malloc_trials_list[k]=atoi(argv[arg_count]);
                }
                else
                {
                    fprintf(stderr, "WARNING: MISSING malloc trial\n");
                    NUM_OF_MALLOC_T=1;
                    malloc_trials_list[0]=1000;//INT_MAX;
                }
            }
        }

        #ifdef SLIP_TCOV_SHOW_LIST
        printf ("gmp ntrials list is: ");
        for (int k=0; k<NUM_OF_TRIALS; k++)
        {
            printf("%ld   ",gmp_ntrial_list[k]);
        }
        printf("\nmalloc trial list is: ");
        for (int k=0; k<NUM_OF_MALLOC_T; k++)
        {
            printf("%d   ",malloc_trials_list[k]);
        }
        printf("\n");
        #endif /* SLIP_TCOV_SHOW_LIST */
    }

    // For SIMPLE_TEST, outter loop iterates for slip_gmp_ntrials initialized
    // from list1 (input for cov_test) and inner loop interates for malloc_count
    // initialized from list2 (input for cov_test)
    //
    // For non SIMPLE_TEST, outter loop iterates for Ab_type from 0 to 5, and
    // set rat correspondingly, and inner loop iterates for malloc_count
    // initialized from 0 to 1000, break when malloc_count>0 at the end of inner
    // loop
    for (int k=0; k<NUM_OF_TRIALS; k++)
    {
        if (IS_SIMPLE_TEST)
        {
            if (k == 1)
            {
                NUM_OF_MALLOC_T=1;
                malloc_trials_list[0]=INT_MAX;
            }
        }
        else
        {
            Ab_type = k;
            //if(k == 1){return 0;} Ab_type = 4;
            rat = rat_list[Ab_type];
            NUM_OF_MALLOC_T = 1000;
        }

        for (int kk=0; kk<NUM_OF_MALLOC_T; kk++)
        {
            if (IS_SIMPLE_TEST)
            {
                slip_gmp_ntrials=gmp_ntrial_list[k];
                printf("initial slip_gmp_ntrials=%ld\n",slip_gmp_ntrials);
                malloc_count=malloc_trials_list[kk];
                printf("%d out of %d, initial malloc_count=%d\n", kk,
                    NUM_OF_MALLOC_T, malloc_count);
            }
            else
            {
                malloc_count = kk;
                printf("[Ab_type rat malloc_count] = [%d %d %d]\n", Ab_type,
                    rat, malloc_count);
            }

            //------------------------------------------------------------------
            // Initialize SLIP LU process
            //------------------------------------------------------------------
            SLIP_initialize();

            //------------------------------------------------------------------
            // Allocate memory
            //------------------------------------------------------------------
            int n=4, numRHS=1, j, nz=11;
            SLIP_info ok;
            SLIP_options* option = SLIP_create_default_options();
            if (!option) {continue;}
            FILE* mat_file = NULL;
            // used in Ab_type = 0
            mpz_t   **B_mpz    = NULL;
            mpz_t   *Ax_mpz    = NULL;
            // used in Ab_type = 1
            double  **B_doub   = NULL;
            double  *Ax_doub   = NULL;
            // used in Ab_type = 2
            int     **B_int    = NULL;
            int     *Ax_int    = NULL;
            // used in Ab_type = 3
            mpq_t   **B_mpq    = NULL;
            mpq_t   *Ax_mpq    = NULL;
            // used in Ab_type = 4
            mpfr_t  **B_mpfr   = NULL;
            mpfr_t  *Ax_mpfr   = NULL;
            // used in rat = 1, 2, 3 correspondingly
            mpq_t   **sol_mpq  = NULL;
            double  **sol_doub = NULL;
            mpfr_t  **sol_mpfr = NULL;
            // used in Ab_type = 5
            mpz_t  *x_mpz  = NULL;
            mpq_t  *x_mpq  = NULL;
            mpfr_t *x_mpfr = NULL;
            double *x_doub = NULL;
            SLIP_sparse *M = NULL;
            // Allocate A
            SLIP_sparse *A = SLIP_create_sparse();
            SLIP_dense *b = SLIP_create_dense();
            // for Column ordering
            SLIP_LU_analysis* S = NULL ;
            option->print_level = 0;

            if (!A || !b) {SLIP_FREE_WORKSPACE; continue;}

            if (Ab_type==0) //mpz
            {
                B_mpz = SLIP_create_mpz_mat(n, numRHS);
                Ax_mpz = SLIP_create_mpz_array(nz);
                if (!B_mpz || !Ax_mpz) {SLIP_FREE_WORKSPACE; continue;}

                for (j = 0; j < n; j++)                           // Get b
                {
                    TEST_CHECK(SLIP_mpz_set_ui(B_mpz[j][0],bxnum3[j]));
                }
                for (j = 0; j < nz; j++)                          // Get Ax
                {
                    TEST_CHECK(SLIP_mpz_set_ui(Ax_mpz[j],Axnum3[j]));
                }

                //failure due to invalid input
                TEST_CHECK_FAILURE(SLIP_build_dense_mpz(b, NULL, n, numRHS));
                TEST_CHECK_FAILURE(SLIP_build_dense_mpz(b, B_mpz, 0, numRHS));
                TEST_CHECK_FAILURE(SLIP_build_sparse_csc_mpz(A, Ap, Ai, NULL,
                    n, nz));

                // successful case
                TEST_CHECK(SLIP_build_sparse_csc_mpz(A, Ap, Ai, Ax_mpz, n, nz));
        TEST_CHECK (SLIP_spok (A, option)) ;
                TEST_CHECK(SLIP_build_dense_mpz(b, B_mpz, n, numRHS));
                option->pivot = SLIP_DIAGONAL;
            }
            else if (Ab_type==1) //double see example4.c
            {
                //failure due to NULL input
                TEST_CHECK_FAILURE(SLIP_build_dense_double(b, NULL, n, numRHS,
            option));
                TEST_CHECK_FAILURE(SLIP_build_sparse_csc_double(A, Ap, Ai, NULL,
                    n, nz, option));

                Ax_doub = (double*) SLIP_calloc(nz, sizeof(double));
                B_doub = SLIP_create_double_mat(n, numRHS);
                if (!B_doub || !Ax_doub) {SLIP_FREE_WORKSPACE; continue;}

                // create empty A and b using uninitialized double mat/array
                TEST_CHECK(SLIP_build_sparse_csc_double(A, Ap, Ai,
                    Ax_doub, n, nz, option));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A
                TEST_CHECK(SLIP_build_dense_double(b, B_doub, n,
                    numRHS, option));
                CLEAR_SLIP_MAT_B;

                // trigger gcd == 1
                for (j = 0; j < n; j++)                           // Get b
                {
                    B_doub[j][0] = bxnum[j]/1e17;
                }
                TEST_CHECK(SLIP_build_dense_double(b, B_doub, n, numRHS,
            option));
                CLEAR_SLIP_MAT_B;

                // trigger gcd != 1
                for (j = 0; j < n; j++)                           // Get b
                {
                    B_doub[j][0] = bxnum[j];
                }
                TEST_CHECK(SLIP_build_dense_double(b, B_doub, n, numRHS,
            option));
                //B_doub[0][0] = 0;
                for (j = 0; j < nz; j++)                          // Get Ax
                {
                    Ax_doub[j] = Axnum[j]/Axden[j];
                }
                TEST_CHECK(SLIP_build_sparse_csc_double(A, Ap, Ai, Ax_doub, n,
                    nz, option));
                TEST_CHECK (SLIP_spok (A, option)) ;
                option->pivot = SLIP_SMALLEST;
            }
            else if (Ab_type==2)//int
            {
                //failure due to NULL input
                TEST_CHECK_FAILURE(SLIP_build_dense_int(b, NULL, n, numRHS));
                TEST_CHECK_FAILURE(SLIP_build_sparse_csc_int(A, Ap, Ai, NULL,
                    n, nz));

                B_int = SLIP_create_int_mat(n, numRHS);
                Ax_int = (int32_t*) SLIP_calloc(nz, sizeof(int32_t));
                if (!B_int || !Ax_int) {SLIP_FREE_WORKSPACE; continue;}

                for (j = 0; j < n; j++)                           // Get b
                {
                    B_int[j][0]=bxnum3[j];
                }
                for (j = 0; j < nz; j++)                          // Get Ax
                {
                    Ax_int[j]=Axnum3[j];
                }

                TEST_CHECK(SLIP_build_sparse_csc_int(A, Ap, Ai, Ax_int, n, nz));
                TEST_CHECK (SLIP_spok (A, option)) ;
                TEST_CHECK(SLIP_build_dense_int(b, B_int, n, numRHS));
            }
            else if (Ab_type==3)//mpq
            {
                //failure due to NULL input
                TEST_CHECK_FAILURE(SLIP_build_dense_mpq(b, NULL, n, numRHS));
                TEST_CHECK_FAILURE(SLIP_build_sparse_csc_mpq(A, Ap, Ai, NULL,
                    n, nz));

                B_mpq = SLIP_create_mpq_mat(n, numRHS);
                Ax_mpq = SLIP_create_mpq_array(nz);
                if (!B_mpq || !Ax_mpq) {SLIP_FREE_WORKSPACE; continue;}

                for (j = 0; j < n; j++)                           // Get b
                {
                    TEST_CHECK(SLIP_mpq_set_ui(B_mpq [j][0], bxnum3[j],
                        bxden3[j]));
                }
                for (j = 0; j < nz; j++)                          // Get Ax
                {
                    TEST_CHECK(SLIP_mpq_set_ui(Ax_mpq [j],Axnum3[j],Axden3[j]));
                }

                TEST_CHECK(SLIP_build_sparse_csc_mpq(A, Ap, Ai, Ax_mpq, n, nz));
                TEST_CHECK (SLIP_spok (A, option)) ;
                TEST_CHECK(SLIP_build_dense_mpq(b, B_mpq, n, numRHS));
            }
            else if (Ab_type==4)//mpfr see example3.c
            {
                //failure due to NULL input
                TEST_CHECK_FAILURE(SLIP_build_dense_mpfr(b, NULL, n, numRHS,
                    option));
                TEST_CHECK_FAILURE(SLIP_build_sparse_csc_mpfr(A, Ap, Ai, NULL,
                    n, nz, option));

                B_mpfr = SLIP_create_mpfr_mat(n, numRHS, option);
                Ax_mpfr = SLIP_create_mpfr_array(nz, option);
                if (!B_mpfr|| !Ax_mpfr) {SLIP_FREE_WORKSPACE; continue;}

                // create empty A and b using uninitialized double mat/array
                TEST_CHECK(SLIP_build_sparse_csc_mpfr(A, Ap, Ai,
                    Ax_mpfr, n, nz, option));
                TEST_CHECK(SLIP_build_dense_mpfr(b, B_mpfr, n,
                    numRHS, option));
                // to trigger SLIP_SINGULAR
                TEST_CHECK(SLIP_LU_analyze(&S, A, option));
                sol_mpq = SLIP_create_mpq_mat(n, numRHS);
                TEST_CHECK_FAILURE(SLIP_solve_mpq(sol_mpq, A, S, b, option));
                option->pivot = SLIP_LARGEST;
                TEST_CHECK_FAILURE(SLIP_solve_mpq(sol_mpq, A, S, b, option));
                option->pivot = SLIP_FIRST_NONZERO;
                TEST_CHECK_FAILURE(SLIP_solve_mpq(sol_mpq, A, S, b, option));
                //free the memory alloc'd
                SLIP_delete_mpq_mat(&sol_mpq, n, numRHS);
                SLIP_delete_LU_analysis(&S) ;
                CLEAR_SLIP_MAT_A;
                CLEAR_SLIP_MAT_B;

                // trigger gcd == 1
                int prec = option->prec;
                option->prec = 17;
                for (j = 0; j < n; j++)                               // Get B
                {
                    TEST_CHECK(SLIP_mpfr_set_d(B_mpfr[j][0], bxnum[j]/1e17,
                        MPFR_RNDN));
                }
                for (j = 0; j < nz; j++)                             // Get Ax
                {
                    TEST_CHECK(SLIP_mpfr_set_d(Ax_mpfr[j], Axnum[j]/1e17,
                        MPFR_RNDN));
                    TEST_CHECK(SLIP_mpfr_div_d(Ax_mpfr[j], Ax_mpfr[j], Axden[j],
                        MPFR_RNDN));
                }
                TEST_CHECK(SLIP_build_dense_mpfr(b, B_mpfr, n, numRHS, option));
                TEST_CHECK(SLIP_build_sparse_csc_mpfr(A, Ap, Ai, Ax_mpfr,
                    n, nz, option));
                CLEAR_SLIP_MAT_A;
                CLEAR_SLIP_MAT_B;
                option->prec = prec;

                // trigger gcd != 1
                for (j = 0; j < n; j++)                               // Get B
                {
                    TEST_CHECK(SLIP_mpfr_set_d(B_mpfr[j][0], bxnum[j],
                        MPFR_RNDN));
                    TEST_CHECK(SLIP_mpfr_div_d(B_mpfr[j][0], B_mpfr[j][0],
                        bxden[j], MPFR_RNDN));
                }
                for (j = 0; j < nz; j++)                             // Get Ax
                {
                    TEST_CHECK(SLIP_mpfr_set_d(Ax_mpfr[j], Axnum[j],MPFR_RNDN));
                    TEST_CHECK(SLIP_mpfr_div_d(Ax_mpfr[j], Ax_mpfr[j], Axden[j],
                        MPFR_RNDN));
                }
                TEST_CHECK(SLIP_build_sparse_csc_mpfr(A, Ap, Ai, Ax_mpfr,
                    n, nz, option));
                TEST_CHECK (SLIP_spok (A, option)) ;
                TEST_CHECK(SLIP_build_dense_mpfr(b, B_mpfr, n, numRHS, option));
            }
            else
            {
                /*
                //test for SLIP_tripread and SLIP_tripread_double
                M = SLIP_create_sparse();

                // fail case
                char *bad_mat1 = "../ExampleMats/bad_mat1.txt";
                mat_file = fopen(bad_mat1,"r");
                if( mat_file == NULL )
                {
                    perror("Error while opening the file");
                    SLIP_FREE_WORKSPACE;
                    continue;
                }
                TEST_CHECK_FAILURE(SLIP_tripread_double(M, mat_file));
                CLEAR_SLIP_MAT_M;                //free the memory alloc'd in M
                fclose (mat_file);
                mat_file = NULL;

                char *bad_mat2 = "../ExampleMats/bad_mat2.txt";
                mat_file = fopen(bad_mat2,"r");
                if( mat_file == NULL )
                {
                    perror("Error while opening the file");
                    SLIP_FREE_WORKSPACE;
                    continue;
                }
                TEST_CHECK_FAILURE(SLIP_tripread_double(M, mat_file));
                CLEAR_SLIP_MAT_M;                //free the memory alloc'd in M
                fclose (mat_file);
                mat_file = NULL;

                // successful case
                char *mat_name = "../ExampleMats/test_mat.txt";
                mat_file = fopen(mat_name,"r");
                if( mat_file == NULL )
                {
                    perror("Error while opening the file");
                    SLIP_FREE_WORKSPACE;
                    continue;
                }
                TEST_CHECK(SLIP_tripread(M, mat_file));
                // move cursor to the beginning of the file
                if ( fseek(mat_file, 0L, SEEK_SET) != 0 )
                {
                    SLIP_FREE_WORKSPACE;
                    continue;
                }
                CLEAR_SLIP_MAT_M;                //free the memory alloc'd in M
                TEST_CHECK(SLIP_tripread_double(M, mat_file));
                fclose (mat_file);
                mat_file = NULL;
                */

                // test for SLIP_build_sparse_trip*
                // and failure of some functions
                n = 4, nz = 11;
                int I[11]={0, 1, 2, 2, 3, 1, 2, 3, 0, 1, 2};
                int J[11]={0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3};

                double x_doub2[11] = {1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4};
                int x_int[11] = {1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4};
                x_mpz  = SLIP_create_mpz_array(nz);
                x_mpq  = SLIP_create_mpq_array(nz);
                x_mpfr = SLIP_create_mpfr_array(nz, option);
                x_doub = (double*) SLIP_calloc(nz, sizeof(double));
                if (!x_mpz || !x_mpq || !x_mpfr || !x_doub)
                {
                    SLIP_FREE_WORKSPACE;
                    continue;
                }

                for (j = 0; j < nz; j++)
                {
                    TEST_CHECK(SLIP_mpq_set_ui(x_mpq[j],2*x_int[j],2));
                    TEST_CHECK(SLIP_mpfr_set_d(x_mpfr[j],x_doub2[j],MPFR_RNDN));
                    TEST_CHECK(SLIP_mpfr_div_d(x_mpfr[j], x_mpfr[j], 1,
                        MPFR_RNDN));
                    x_doub[j] = x_doub2[j];
                    TEST_CHECK(SLIP_mpz_set_si(x_mpz[j], x_int[j]));
                }

                // failure case: incorrect index array input
                int I1[11]={0, -1, -2, 2, 3, 1, 2, 3, 0, 1, 2};
                int P1[11]={0, 3, 5, 8, 11};
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpz(A, I1, J, x_mpz,
                    n, nz));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A
                TEST_CHECK_FAILURE(SLIP_build_sparse_csc_mpz(A, P1, I1, x_mpz,
                    n, nz));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A

                // failure cases: input NULL pointer
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpz   (A, I, J, NULL,
                    n, nz));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_double(A, I, J, NULL,
                    n, nz, option));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_int   (A, I, J, NULL,
                    n, nz));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpq   (A, I, J, NULL,
                    n, nz));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpfr  (A, I, J, NULL,
                    n, nz, option));

                // successful cases
                TEST_CHECK(SLIP_build_sparse_trip_mpz(A, I, J, x_mpz, n, nz));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A
                TEST_CHECK(SLIP_build_sparse_trip_double(A, I, J, x_doub, n,
                    nz, option));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A
                TEST_CHECK(SLIP_build_sparse_trip_int(A, I, J, x_int, n, nz));
                ok = SLIP_gmp_printf("scale = %Qd\n",A->scale);
                if (ok < 0) {TEST_CHECK(ok);}
                option->print_level = 3;
                TEST_CHECK(SLIP_spok (A, option));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A
                TEST_CHECK(SLIP_build_sparse_trip_mpq(A, I, J, x_mpq, n, nz));
                ok = SLIP_gmp_printf("scale = %Qd\n",A->scale);
                if (ok < 0) {TEST_CHECK(ok);}
                TEST_CHECK(SLIP_spok (A, option));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A
                TEST_CHECK(SLIP_build_sparse_trip_mpfr(A, I, J, x_mpfr, n, nz,
                    option));
                ok = SLIP_gmp_printf("scale = %Qd\n",A->scale);
                if (ok < 0) {TEST_CHECK(ok);}
                TEST_CHECK(SLIP_spok (A, option));

                //test coverage for SLIP_spok()
                A->i[0] = -1;
                option->print_level = 1;
                TEST_CHECK_FAILURE(SLIP_spok(A, option));
                SLIP_delete_mpz_array(&(A->x), A->nzmax);
                TEST_CHECK_FAILURE(SLIP_spok(A, option));
                A->p[1] = 2;
                A->p[2] = 1;
                TEST_CHECK_FAILURE(SLIP_spok(A, option));
                A->p[0] = 1;
                TEST_CHECK_FAILURE(SLIP_spok(A, option));
                A->nzmax = -1;
                TEST_CHECK_FAILURE(SLIP_spok(A, option));
                A->n = -1;
                TEST_CHECK_FAILURE(SLIP_spok(A, option));
                A->m = -1;
                TEST_CHECK_FAILURE(SLIP_spok(A, option));
                CLEAR_SLIP_MAT_A;                //free the memory alloc'd in A


                //test coverage for slip_gmp_reallocate()
                void *p_new = NULL;
                TEST_CHECK(slip_gmp_realloc_test(&p_new, NULL , 0, 1));
                TEST_CHECK(slip_gmp_realloc_test(&p_new, p_new, 1, 0));
                printf("test\n");

                // Incorrect calling with NULL pointer(s)
                TEST_CHECK_FAILURE(SLIP_LU_analyze(NULL, NULL, NULL));
                TEST_CHECK_FAILURE(SLIP_LU_analyze(&S, NULL, NULL));
                TEST_CHECK_FAILURE(SLIP_LU_factorize(NULL, NULL, NULL, NULL,
                    NULL, NULL, NULL));
                TEST_CHECK_FAILURE(SLIP_LU_solve(NULL, NULL, NULL, NULL, NULL,
                    NULL));

                SLIP_FREE_WORKSPACE;

                // for miscellaneous test, continue to next loop directly
                if (!IS_SIMPLE_TEST)
                {
                    if (malloc_count > 0)
                    {
                        rat_list[5] = kk;
                        break;
                    }
                    else {continue;}
                }
                else
                {
                    continue;
                }
            }

            // Column ordering using either AMD, COLAMD or nothing
            TEST_CHECK(SLIP_LU_analyze(&S, A, option));
            option->print_level = 3;
            int check2;

            //------------------------------------------------------------------
            // SLIP LU Factorization, Solve and verification
            //------------------------------------------------------------------
            if (rat==1)
            {
                sol_mpq = SLIP_create_mpq_mat(n, numRHS);
                TEST_CHECK(SLIP_solve_mpq(sol_mpq, A, S, b, option));
                if (Ab_type == 0)
                {
                    //failure case: NULL pointer input
                    TEST_CHECK_FAILURE(SLIP_check_solution(A, NULL, b));
                    TEST_CHECK_FAILURE(SLIP_get_double_soln(NULL, sol_mpq, n,
                        numRHS));
                    TEST_CHECK_FAILURE(SLIP_get_mpfr_soln(NULL, sol_mpq, n,
                        numRHS, option));

                    TEST_CHECK(SLIP_check_solution(A, sol_mpq, b));
                    check2 = ok;  // track the status of SLIP_check_solution
                    TEST_CHECK(SLIP_print_stats_mpq(stdout, sol_mpq, n, numRHS,
                        check2,option));
                    option->print_level = 3;
                    TEST_CHECK(SLIP_spok (A, option));
                    //SLIP_PRINT_OK(ok);

                    //intentionally change b to fail SLIP_LU_Check()
                    TEST_CHECK(SLIP_mpz_set_ui (b->x[0][0], 1000));
                    check2=SLIP_check_solution(A, sol_mpq, b);

                    // Print result using SLIP_print_stats, which should return
                    // SLIP_INCORRECT since check2 == SLIP_INCORRECT
                    ok=SLIP_print_stats_mpq(stdout, sol_mpq, n, numRHS, check2,
                        option);
                    if (ok!=SLIP_INCORRECT)
                    {
                        SLIP_PRINT_OK(ok);
                        SLIP_FREE_WORKSPACE;
                        continue;
                    }
                }
            }
            else if (rat==2)
            {
                sol_doub = SLIP_create_double_mat(n, numRHS);
                TEST_CHECK(SLIP_solve_double(sol_doub, A, S, b, option));
            }
            else
            {
                sol_mpfr = SLIP_create_mpfr_mat(n, numRHS, option);
                TEST_CHECK(SLIP_solve_mpfr(sol_mpfr, A, S, b, option));
                TEST_CHECK(SLIP_print_stats_mpfr(stdout, sol_mpfr, n, numRHS,
                    SLIP_OK, option));
            }
            //------------------------------------------------------------------
            // Free Memory
            //------------------------------------------------------------------
            SLIP_FREE_WORKSPACE;
            if(!IS_SIMPLE_TEST)
            {
                if (malloc_count > 0)
                {
                    rat_list[k] = kk;
                    break;
                }
                else {continue;}
            }
        }
    }
    if (IS_SIMPLE_TEST)
    {
        free(gmp_ntrial_list);
        free(malloc_trials_list);
    }
    else
    {
        printf("malloc_count for Ab_type = 0 ~ 5 are %d %d %d %d %d %d\n",
            rat_list[0], rat_list[1], rat_list[2], rat_list[3], rat_list[4],
            rat_list[5]);
    }
    return 0;
}

