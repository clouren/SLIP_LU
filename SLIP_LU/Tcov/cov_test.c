//------------------------------------------------------------------------------
// SLIP_LU/Tcov/cov_test.c: test coverage for SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/*
 * When the test is run without input argument, brutal test is used and simple
 * test otherwise. Read the following for detailed instruction and information
 *
 * For simple test, the test needs to be run with command
 * ./cov_test Ab_type xtype N list1[0] ... list1[N-1] M list2[0] ... list2[M-1]
 * Ab_type: type of Matrix A and vector b: 0 mpz, 1 double, 2 int64, 3 mpq,
 *      4 mpfr, 5 for miscellaneous test.
 * xtype: type of solution: 1: full precision rational arithmetic,
 *                        2: double, 3:  user specified precision.
 * N and list1 specify the test list for slip_gmp_ntrials (in SLIP_gmp.h)
 * M and list2 specify the test list for malloc_count (in tcov_malloc_test.h)
 * N, list1, M, list2 are optional, but N and list1 are required when M and
 * list2 is wanted
 *
 * For brutal test, the test is run with command
 * ./cov_test
 * the test will run through all cases
 * (specifically, [Ab_type xtype]={[0 1], [1 1], [2 2], [3 3], [4 3], [5 3]})
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

#define SLIP_FREE_ALL                            \
{                                                \
    SLIP_LU_analysis_free(&S);                   \
    SLIP_matrix_free(&A,option);                 \
    SLIP_matrix_free(&b, option);                \
    SLIP_matrix_free(&B, option);                \
    SLIP_matrix_free(&Ax, option);               \
    SLIP_matrix_free(&sol, option);              \
    SLIP_FREE(option);                           \
    if (mat_file != NULL) {fclose(mat_file);}    \
    SLIP_finalize() ;                            \
}

#include "SLIP_LU_internal.h"

#define TEST_CHECK(method)                       \
{                                                \
    info = (method) ;                            \
    if (info != SLIP_OK)                         \
    {                                            \
        SLIP_PRINT_INFO (info) ;                 \
        SLIP_FREE_ALL;                     \
        continue;                                \
    }                                            \
}

#define TEST_CHECK_FAILURE(method)               \
{                                                \
    info = (method) ;                            \
    if (info != SLIP_INCORRECT_INPUT && info != SLIP_SINGULAR) \
    {                                            \
        SLIP_PRINT_INFO (info) ;                 \
        SLIP_FREE_ALL ;                    \
        continue ;                               \
    }                                            \
    else                                         \
    {                                            \
        printf("Expected failure at line %d\n", __LINE__);\
    }                                            \
}

#define MAX_MALLOC_COUNT 10000

int64_t Ap[5] = {0, 3, 5, 8, 11};
int64_t Ai[11]   = {0, 1, 2, 2, 3, 1, 2, 3, 0, 1,  2};
double Axnum[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};  // Numerator of x
double Axden[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};  // Denominator of x
double bxnum[4] = {170, 1820, 61, 670};                // Numerator of b
double bxden[4] = {15,  3,   6,  7};                    // Denominator of b
int64_t Axnum3[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};    // Numerator of x
int64_t Axden3[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};    // Denominator of x
int64_t bxnum3[4] = {17, 182, 61, 67};                      // Numerator of b
int64_t bxden3[4] = {15,  3,   6,  7};                      // Denominator of b

int main( int argc, char* argv[])
{
    bool IS_SIMPLE_TEST = true;
    int rat_list[6] = {1, 1, 2, 3, 3, 3};  // only used in brutal test
    int Ab_type = 0, xtype = 1 ;
    int64_t NUM_OF_TRIALS = 0 ;
    int64_t NUM_OF_MALLOC_T = 0;
    int64_t *gmp_ntrial_list=NULL;         // only used in simple test
    int64_t *malloc_trials_list=NULL;          // only used in simple test

    //--------------------------------------------------------------------------
    // parse input arguments
    //--------------------------------------------------------------------------

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

        int64_t arg_count = 0;
        // type of Matrix A and vector b:
        // 0 mpz, 1 double, 2 int64_t, 3 mpq, 4 mpfr
        Ab_type = atoi(argv[++arg_count]);
        //type of solution: 1: full precision rational arithmetic,
        //                  2: double, 3:  user specified precision.
        xtype=atoi(argv[++arg_count]);
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
            for (int64_t k=0; k<NUM_OF_TRIALS; k++)
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
            malloc_trials_list=SLIP_malloc(NUM_OF_MALLOC_T* sizeof(int64_t));
            malloc_trials_list[0]=MAX_MALLOC_COUNT;//INT_MAX;
        }
        else
        {
            NUM_OF_MALLOC_T=atoi(argv[arg_count]);
            malloc_trials_list=SLIP_malloc(NUM_OF_MALLOC_T* sizeof(int64_t));
            for (int64_t k=0; k<NUM_OF_MALLOC_T; k++)
            {
                if (argv[++arg_count])
                {
                    malloc_trials_list[k]=atoi(argv[arg_count]);
                }
                else
                {
                    fprintf(stderr, "WARNING: MISSING malloc trial\n");
                    NUM_OF_MALLOC_T=1;
                    malloc_trials_list[0]=MAX_MALLOC_COUNT;//INT_MAX;
                }
            }
        }

        #ifdef SLIP_TCOV_SHOW_LIST
        printf ("gmp ntrials list is: ");
        for (int64_t k=0; k<NUM_OF_TRIALS; k++)
        {
            printf("%ld   ",gmp_ntrial_list[k]);
        }
        printf("\nmalloc trial list is: ");
        for (int64_t k=0; k<NUM_OF_MALLOC_T; k++)
        {
            printf("%d   ",malloc_trials_list[k]);
        }
        printf("\n");
        #endif /* SLIP_TCOV_SHOW_LIST */
    }

    //--------------------------------------------------------------------------
    // run all trials
    //--------------------------------------------------------------------------

    // For SIMPLE_TEST, outter loop iterates for slip_gmp_ntrials initialized
    // from list1 (input for cov_test) and inner loop interates for
    // malloc_count initialized from list2 (input for cov_test).
    //
    // For non SIMPLE_TEST, outter loop iterates for Ab_type from 0 to 5, and
    // set xtype correspondingly, and inner loop iterates for malloc_count
    // initialized from 0 to MAX_MALLOC_COUNT, break when malloc_count>0 at the
    // end of inner loop.

    for (int64_t k=0; k<NUM_OF_TRIALS; k++)
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
            xtype = rat_list[Ab_type];
            NUM_OF_MALLOC_T = MAX_MALLOC_COUNT;
        }

        for (int64_t kk=0; kk<NUM_OF_MALLOC_T; kk++)
        {
            if (IS_SIMPLE_TEST)
            {
                slip_gmp_ntrials=gmp_ntrial_list[k];
                printf("initial slip_gmp_ntrials=%ld\n",slip_gmp_ntrials);
                malloc_count=malloc_trials_list[kk];
                printf("%"PRId64" out of %"PRId64", "
                    "initial malloc_count=%"PRId64"\n",
                    kk, NUM_OF_MALLOC_T, malloc_count);
            }
            else
            {
                malloc_count = kk;
                printf("[Ab_type xtype malloc_count] = "
                    "[%d %d %"PRId64"]\n", Ab_type, xtype, malloc_count);
            }

            //------------------------------------------------------------------
            // Initialize SLIP LU process
            //------------------------------------------------------------------

            SLIP_initialize();

            //------------------------------------------------------------------
            // Allocate memory
            //------------------------------------------------------------------

            int64_t n=4, numRHS=1, j, nz=11;
            SLIP_info info ;
            SLIP_options* option = SLIP_create_default_options();
            if (!option) {continue;}
            option->print_level = 3;
            FILE* mat_file = NULL;

            // used in as source in different Ab_type for A and b
            SLIP_matrix *B   = NULL;
            SLIP_matrix *Ax  = NULL;

            // matrix A, b and solution
            SLIP_matrix *A = NULL ;
            SLIP_matrix *b = NULL ;
            SLIP_matrix *sol = NULL;

            // for Column ordering
            SLIP_LU_analysis* S = NULL ;

            if (Ab_type >= 0 && Ab_type <= 4)
            {

                //--------------------------------------------------------------
                // Solve A*x=b where A and b are created from mpz entries
                //--------------------------------------------------------------

                TEST_CHECK(SLIP_matrix_allocate(&B, SLIP_DENSE,
                    (SLIP_type) Ab_type, n,
                    numRHS, n*numRHS, false, true, option));
                TEST_CHECK(SLIP_matrix_allocate(&Ax, SLIP_CSC,
                    (SLIP_type) Ab_type, n,
                    n, nz, false, true, option));

                // fill Ax->i and Ax->p
                for (j = 0; j < n+1; j++)
                {
                    Ax->p[j] = Ap[j];
                }
                for (j = 0; j < nz; j++)
                {
                    Ax->i[j] = Ai[j];
                }

                // special failure cases
                if (Ab_type == 2)// MPFR
                {
                    // create empty A and b using uninitialized double mat/array
                    // to trigger all-zero array condition
                    TEST_CHECK(SLIP_matrix_copy(&A, SLIP_CSC, SLIP_MPZ, Ax,
                        option));
                    TEST_CHECK(SLIP_matrix_copy(&b, SLIP_DENSE, SLIP_MPZ, B,
                        option));
                    // to trigger SLIP_SINGULAR
                    TEST_CHECK_FAILURE(SLIP_backslash(&sol, SLIP_MPQ, A, b,
                        option));
                    option->pivot = SLIP_LARGEST;
                    TEST_CHECK_FAILURE(SLIP_backslash(&sol, SLIP_MPQ, A, b,
                        option));
                    option->pivot = SLIP_FIRST_NONZERO;
                    TEST_CHECK_FAILURE(SLIP_backslash(&sol, SLIP_MPQ, A, b,
                       option));
                    // incorrect solution type
                    TEST_CHECK_FAILURE(SLIP_backslash(&sol, SLIP_MPZ, A, b,
                       option));

                    // A->p = NULL
                    SLIP_FREE(A->p);
                    TEST_CHECK_FAILURE(SLIP_backslash(&sol, SLIP_MPQ, A, b,
                       option));

                    //free the memory alloc'd
                    SLIP_matrix_free (&A, option) ;
                    SLIP_matrix_free (&b, option) ;

                    // trigger gcd == 1
                    uint64_t prec = option->prec;
                    option->prec = 17;
                    double pow2_17 = pow(2,17);
                    for (j = 0; j < n; j++)                             // Get B
                    {
                        TEST_CHECK(SLIP_mpfr_set_d(SLIP_2D(B,j,0,mpfr),
                            bxnum[j]/pow2_17, MPFR_RNDN));
                    }
                    TEST_CHECK (SLIP_matrix_check (B, option)) ;
                    TEST_CHECK(SLIP_matrix_copy(&b, SLIP_DENSE, SLIP_MPZ,B,
                        option));
                    SLIP_matrix_free (&A, option) ;
                    SLIP_matrix_free (&b, option) ;

                    // restore default precision
                    option->prec = prec;
                }
                else if (Ab_type == 4)// double
                {
                    // create empty A using uninitialized double mat/array
                    // to trigger all-zero array condition
                    TEST_CHECK(SLIP_matrix_copy(&A, SLIP_CSC, SLIP_MPZ, Ax,
                        option));

                    // trigger gcd == 1
                    for (j = 0; j < n; j++)                           // Get b
                    {
                        SLIP_2D(B,j,0,fp64) = bxnum[j]/1e17;
                    }
                    TEST_CHECK(SLIP_matrix_copy(&b, SLIP_DENSE, SLIP_MPZ, B,
                        option));
                    SLIP_matrix_free (&b, option) ;

                    // failure case: Ax->x = NULL
                    double *tmp_Ax_fp64 = Ax->x.fp64;
                    Ax->x.fp64 = NULL;
                    TEST_CHECK_FAILURE(SLIP_matrix_copy(&A, SLIP_CSC, SLIP_MPZ,
                        Ax,option));
                    Ax->x.fp64 = tmp_Ax_fp64;
                }

                // fill Ax->x and b->x
                for (j = 0; j < n; j++)                           // Get b
                {
                    if (Ab_type == 0) //MPZ
                    {
                        TEST_CHECK(SLIP_mpz_set_ui(SLIP_2D(B, j, 0, mpz),
                            bxnum3[j]));
                    }
                    else if (Ab_type == 1)// MPQ
                    {
                        TEST_CHECK(SLIP_mpq_set_ui(SLIP_2D(B,j,0,mpq),
                            bxnum3[j], bxden3[j]));
                    }
                    else if (Ab_type == 2)// MPFR
                    {
                        TEST_CHECK(SLIP_mpfr_set_d(SLIP_2D(B,j,0,mpfr),bxnum[j],
                            MPFR_RNDN));
                        TEST_CHECK(SLIP_mpfr_div_d(SLIP_2D(B,j,0,mpfr),
                            SLIP_2D(B,j,0,mpfr), bxden[j], MPFR_RNDN));
                    }
                    else if (Ab_type == 3)// INT64
                    {
                        SLIP_2D(B,j,0,int64)=bxnum3[j];
                    }
                    else // double
                    {
                        SLIP_2D(B,j,0,fp64) = bxnum[j];
                    }
                }
                for (j = 0; j < nz; j++)                          // Get Ax
                {
                    if (Ab_type == 0)
                    {
                        TEST_CHECK(SLIP_mpz_set_ui(Ax->x.mpz[j],Axnum3[j]));
                    }
                    else if (Ab_type == 1)
                    {
                        TEST_CHECK(SLIP_mpq_set_ui(Ax->x.mpq[j],Axnum3[j],
                            Axden3[j]));
                    }
                    else if (Ab_type == 2)
                    {
                        TEST_CHECK(SLIP_mpfr_set_d(Ax->x.mpfr[j], Axnum[j],
                            MPFR_RNDN));
                        TEST_CHECK(SLIP_mpfr_div_d(Ax->x.mpfr[j], Ax->x.mpfr[j],
                            Axden[j], MPFR_RNDN))
                    }
                    else if (Ab_type == 3)
                    {
                        Ax->x.int64[j]=Axnum3[j];
                    }
                    else
                    {
                        Ax->x.fp64[j] = Axnum[j]/Axden[j];
                    }
                }

                // successful case
                TEST_CHECK(SLIP_matrix_copy(&A, SLIP_CSC, SLIP_MPZ, Ax,option));
                TEST_CHECK(SLIP_matrix_copy(&b, SLIP_DENSE, SLIP_MPZ,B,option));
                TEST_CHECK (SLIP_matrix_check (A, option)) ;
                TEST_CHECK (SLIP_matrix_check (b, option)) ;
                if (Ab_type == 0)
                {
                    option->pivot = SLIP_DIAGONAL;
                }
                else if (Ab_type == 1)
                {
                    option->pivot = SLIP_SMALLEST;
                }
            }
            else
            {

                //--------------------------------------------------------------
                // create A and b from triplets
                //--------------------------------------------------------------

                // test for SLIP_build_sparse_trip*
                // and failure of some functions
                n = 4, nz = 11;
                int64_t m1, n1, nz1;
                int64_t I[11]={0, 1, 2, 2, 3, 1, 2, 3, 0, 1, 2};
                int64_t J[11]={0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3};
                int64_t P[11]={0, 3, 5, 8, 11};
                // failure case: incorrect index array input
                int64_t I1[11]={0, -1, -2, 2, 3, 1, 2, 3, 0, 1, 2};

                double x_doub2[11] = {1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4};
                int64_t x_int64[11] = {1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4};
                for (int tk = 0; tk < 15; tk++)
                {
                    int type = tk%5;
                    int kind = tk/5;
                    if (kind != 2)
                    {
                        m1 = n;
                        n1 = n;
                        nz1 = nz;
                    }
                    else
                    {
                        m1 = nz1;
                        n1 = 1;
                        nz1 = nz1;
                    }
                    TEST_CHECK(SLIP_matrix_allocate(&Ax, (SLIP_type) kind,
                        (SLIP_type)type, m1, n1, nz1, false, true, option));

                    // fill Ax->p
                    if(kind == 0)
                    {
                        for (j = 0; j < n+1; j++)
                        {
                            Ax->p[j] = P[j];
                        }
                    }
                    // fill Ax->i and Ax->j
                    for (j = 0; j < nz; j++)
                    {
                        if (kind != 2) {Ax->i[j] = I[j];}
                        // triplet
                        if (kind == 1){  Ax->j[j] = J[j];}
                        switch (type)
                        {
                            case 0: // MPZ
                            {
                                TEST_CHECK(SLIP_mpz_set_si(Ax->x.mpz[j],
                                    x_int64[j]));
                            }
                            break;

                            case 1: // MPQ
                            {
                                TEST_CHECK(SLIP_mpq_set_ui(Ax->x.mpq[j],
                                    2*x_int64[j],2));
                            }
                            break;

                            case 2: // MPFR
                            {
                                TEST_CHECK(SLIP_mpfr_set_d(Ax->x.mpfr[j],
                                    x_doub2[j], MPFR_RNDN));
                                TEST_CHECK(SLIP_mpfr_div_d(Ax->x.mpfr[j],
                                    Ax->x.mpfr[j], 1, MPFR_RNDN));
                            }
                            break;

                            case 3: // INT64
                            {
                                Ax->x.int64[j] = x_int64[j];
                            }
                            break;

                            case 4: // double
                            {
                                Ax->x.fp64[j] = x_doub2[j];
                            }
                            break;

                            default: break;
                        }
                    }
                    TEST_CHECK (SLIP_matrix_check (Ax, option));

                    // convert to all different type of matrix
                    for (int tk1 = 0; tk1 < 15; tk1++)
                    {
                        // successful cases
                        int type1 = tk1%5;
                        int kind1 = tk1/5;
                        printf("converting from %s(%d) %s(%d) to %s(%d) "
                            "%s(%d)\n",kind < 1 ? "CSC" : kind < 2 ? "Triplet" :
                            "Dense",kind, type < 1 ? "MPZ" : type < 2 ? "MPQ" :
                            type <3 ?  "MPFR" : type < 4 ? "int64" :
                            "double",type, kind1 < 1 ?  "CSC" : kind1 < 2 ?
                            "Triplet" : "Dense", kind1,type1 < 1 ? "MPZ" :
                            type1 < 2 ? "MPQ" : type1 < 3 ?  "MPFR" : type1 < 4
                            ? "int64" : "double",type1) ;

                        TEST_CHECK(SLIP_matrix_copy(&A, (SLIP_kind)kind1,
                            (SLIP_type) type1, Ax, option));

                        // just perform once to perform some failure cases
                        if (tk == 0 && tk1 == 0)
                        {
                            // test SLIP_matrix_check
                            A->i[0] = -1;
                            TEST_CHECK_FAILURE(SLIP_matrix_check(A, option));
                            for (int64_t i = 0; i < A->nzmax; i++) 
                            { 
                                if ( A->x.mpz[i] != NULL) 
                                { 
                                    SLIP_MPZ_CLEAR( A->x.mpz[i]); 
                                } 
                            }

                            TEST_CHECK_FAILURE(SLIP_matrix_check(A, option));
                            A->p[1] = 2;
                            A->p[2] = 1;
                            TEST_CHECK_FAILURE(SLIP_matrix_check(A, option));
                            A->p[0] = 1;
                            TEST_CHECK_FAILURE(SLIP_matrix_check(A, option));
                            A->nzmax = -1;
                            TEST_CHECK_FAILURE(SLIP_matrix_check(A, option));
                            A->n = -1;
                            TEST_CHECK_FAILURE(SLIP_matrix_check(A, option));
                            A->m = -1;
                            TEST_CHECK_FAILURE(SLIP_matrix_check(A, option));

                            //test coverage for slip_gmp_reallocate()
                            void *p_new = NULL;
                            TEST_CHECK(slip_gmp_realloc_test(&p_new, NULL,0,1));
                            TEST_CHECK(slip_gmp_realloc_test(&p_new,p_new,1,0));
                            printf("test\n");

                            // Incorrect calling with NULL pointer(s)
                            TEST_CHECK_FAILURE(SLIP_LU_analyze(NULL,NULL,NULL));
                            TEST_CHECK_FAILURE(SLIP_LU_analyze(&S, NULL, NULL));
                            TEST_CHECK_FAILURE(SLIP_LU_factorize(NULL, NULL,
                                NULL, NULL, A, NULL, NULL));
                            SLIP_matrix *L, *U, *rhos;
                            int64_t *pinv ;
                            TEST_CHECK_FAILURE(SLIP_LU_factorize(&L, &U, &rhos,
                                &pinv, A, NULL, NULL));

                            TEST_CHECK(SLIP_matrix_allocate(&b, SLIP_DENSE,
                                SLIP_MPZ, 1, 1, 1, true, true, option));
                            TEST_CHECK(SLIP_matrix_allocate(&L, SLIP_CSC,
                                SLIP_MPZ, 1, 1, 1, true, true, option));
                            TEST_CHECK(SLIP_matrix_allocate(&U, SLIP_CSC,
                                SLIP_MPZ, 1, 1, 1, true, true, option));
                            TEST_CHECK(SLIP_matrix_allocate(&rhos, SLIP_DENSE,
                                SLIP_MPZ, 1, 1, 1, true, true, option));
                            TEST_CHECK_FAILURE(SLIP_LU_solve(NULL, b, A, L, U,
                                rhos, NULL, pinv, option));
                            SLIP_matrix_free (&b, option) ;
                            SLIP_matrix_free (&L, option) ;
                            SLIP_matrix_free (&U, option) ;
                            SLIP_matrix_free (&rhos, option) ;
                            TEST_CHECK_FAILURE(SLIP_LU_solve(NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL));
                        }
                        SLIP_matrix_free (&A, option) ;
                    }
                    SLIP_matrix_free (&Ax, option) ;
                }

#if 0
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpz(&A, I1, J, x_mpz,
                    n, nz));
                SLIP_delete_sparse (&A) ;

                TEST_CHECK_FAILURE(SLIP_build_sparse_csc_mpz(&A, P1, I1, x_mpz,
                    n, nz));
                SLIP_delete_sparse (&A) ;

                // failure cases: input NULL pointer
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpz (&A, I, J, NULL,
                    n, nz));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_double(&A, I, J, NULL,
                    n, nz, option));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_int64 (&A, I, J, NULL,
                    n, nz));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpq   (&A, I, J, NULL,
                    n, nz));
                TEST_CHECK_FAILURE(SLIP_build_sparse_trip_mpfr  (&A, I, J, NULL,
                    n, nz, option));
#endif


                SLIP_FREE_ALL;

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

            //------------------------------------------------------------------
            // SLIP LU backslash
            // solve Ax=b in full precision rational arithmetic
            //------------------------------------------------------------------

            GOTCHA;
            TEST_CHECK(SLIP_backslash(&sol, SLIP_MPQ, A, b, option));
            GOTCHA;

            //------------------------------------------------------------------
            // SLIP LU Factorization, Solve and verification
            //------------------------------------------------------------------

            if (xtype==1)
            {

                if (Ab_type == 4)
                {
                    // This would return SLIP_INCORRECT since sol has been
                    // scaled down so that sol->scale = 1. Therefore sol is
                    // solution for original unscaled Ax=b, while this is
                    // checking if x is the solution for scaled Ax=b
                    TEST_CHECK(slip_check_solution(A, sol, b, option));
                }

            }
            else if (xtype==2)
            {

                //--------------------------------------------------------------
                // copy sol to double
                //--------------------------------------------------------------

                SLIP_matrix *sol_doub;
                TEST_CHECK(SLIP_matrix_copy(&sol_doub, SLIP_DENSE, SLIP_FP64,
                    sol, option));
                SLIP_matrix_free(&sol_doub, option);

            }
            else
            {

                //--------------------------------------------------------------
                // copy sol to mpfr with user-specified precision
                //--------------------------------------------------------------

                SLIP_matrix *sol_mpfr;
                TEST_CHECK(SLIP_matrix_copy(&sol_mpfr, SLIP_DENSE, SLIP_MPFR,
                    sol, option));
                SLIP_matrix_free(&sol_mpfr, option);
            }

            //------------------------------------------------------------------
            // Free Memory
            //------------------------------------------------------------------

            SLIP_FREE_ALL;
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

    //--------------------------------------------------------------------------
    // wrapup
    //--------------------------------------------------------------------------

    if (IS_SIMPLE_TEST)
    {
        free(gmp_ntrial_list);
        free(malloc_trials_list);
    }
    else
    {
        printf("remaining malloc_count for Ab_type = 0~5 are %d %d %d %d "
            "%d %d\n", rat_list[0], rat_list[1], rat_list[2], rat_list[3],
            rat_list[4], rat_list[5]);
    }
    return 0;
}

