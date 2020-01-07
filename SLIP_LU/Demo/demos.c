# include "demos.h"

/* Purpose: This function prints out the user specified/default options
 * this is primarily intended for debugging
 */
void SLIP_print_options // display specified/default options to user
(
    SLIP_options* option // struct containing all of the options
)
{
    char *piv, *order;
    if (option->order == SLIP_COLAMD)
    {
        order = "the COLAMD";
    }
    else if (option->order == SLIP_AMD)
    {
        order = "the AMD";
    }
    else if (option->order == SLIP_NO_ORDERING)
    {
        order = "No";
    }
    else
    {
        order = "(undefined)";
    }

    if (option->pivot == SLIP_SMALLEST)
    {
        piv = "smallest";
    }
    else if (option->pivot == SLIP_DIAGONAL)
    {
        piv = "diagonal";
    }
    else if (option->pivot == SLIP_FIRST_NONZERO)
    {
        piv = "first nonzero";
    }
    else if (option->pivot == SLIP_TOL_SMALLEST)
    {
        piv = "diagonal small tolerance";
    }
    else if (option->pivot == SLIP_TOL_LARGEST)
    {
        piv = "diagonal large tolerance";
    }
    else
    {
        piv = "largest";
    }

    printf("\n\n****COMMAND PARAMETERS****");
    printf("\nUsing %s ordering and selecting the %s pivot", order, piv);
    if (option->pivot == SLIP_TOL_SMALLEST ||
        option->pivot == SLIP_TOL_LARGEST)
    {
        printf("\nTolerance used: %lf\n",option->tol);
    }
}


/* Purpose: This processes the command line for user specified options */
SLIP_info SLIP_process_command_line //processes the command line
(
    int32_t argc,           // number of command line arguments
    char* argv[],           // set of command line arguments
    SLIP_options* option,   // struct containing the command options
    char** mat_name,        // Name of the matrix to be read in
    char** rhs_name,        // Name of the RHS vector to be read in
    int32_t *rat            // data type of output solution.
                            // 1: mpz, 2: double, 3: mpfr
)
{
    for (int32_t i = 1; i < argc; i++)
    {
        char* arg = (char*) argv[i];
        if ( strcmp(arg,"help") == 0)
        {
            SLIP_show_usage();
            return SLIP_INCORRECT_INPUT;
        }
        else if ( strcmp(arg,"p") == 0 || strcmp(arg,"piv") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be a pivot argument between"
                    " 0-5 following p\n");
                return SLIP_INCORRECT_INPUT;
            }
            option->pivot = atoi(argv[i]);
            if (option->pivot < 0 || option->pivot > 5)
            {
                printf("\n****ERROR! Invalid pivot selection!"
                    "\nDefaulting to smallest pivot\n\n");
                option->pivot = SLIP_SMALLEST;
            }
        }
        else if ( strcmp(arg, "q") == 0 || strcmp(arg,"col") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be an argument between 0-2"
                    "following q\n");
                return SLIP_INCORRECT_INPUT;
            }
            option->order = atoi(argv[i]);
            if (option->order < 0 || option->order > 2)
            {
                printf("\n****ERROR! Invalid column ordering"
                    "\nDefaulting to COLAMD\n\n");
                option->order = SLIP_COLAMD;
            }
        }
        else if ( strcmp(arg,"t") == 0 || strcmp(arg, "tol") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be a non negative tolerance"
                    " value following t\n");
                return SLIP_INCORRECT_INPUT;
            }
            else if (!atof(argv[i]))
            {
                printf("\n****ERROR! There must be a non negative tolerance"
                    " value following t\n");
                return SLIP_INCORRECT_INPUT;
            }
            option->tol = atof(argv[i]);
            if (option->tol < 0)
            {
                printf("\n****ERROR! Invalid Tolerance, tolerance must be"
                    " non-negative\n");
                return SLIP_INCORRECT_INPUT;
            }
        }
        else if ( strcmp(arg,"out2") == 0 || strcmp(arg, "o2") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (print nothing) 1 (print err) or 2 (terse) \n");
                return SLIP_INCORRECT_INPUT;
            }
            else if (!atoi(argv[i]))
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (print nothing) 1 (print err) or 2 (terse) \n");
                return SLIP_INCORRECT_INPUT;
            }
            option->print_level = atoi(argv[i]);
        }
        else if ( strcmp(arg, "out") == 0 || strcmp(arg, "o") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! o or out must be followed by"
                    " 1 (rational) 2 (double) or 3 (variable precision) \n");
                return SLIP_INCORRECT_INPUT;
            }
            else if (!atoi(argv[i]))
            {
                printf("\n****ERROR! o or out must be followed by"
                    " 1 (rational) 2 (double) or 3 (variable precision) \n");
                return SLIP_INCORRECT_INPUT;
            }
            *rat = atoi(argv[i]);
            if (*rat < 1 || *rat > 3)
            {
                printf("\n\n****ERROR! Invalid output type!\n"
                   "Defaulting to rational");
                *rat = 1;
            }
            if (*rat == 3)
            {
                if (!argv[++i])
                {
                    printf("\n****ERROR! Precision must be specified\n");
                    return SLIP_INCORRECT_INPUT;
                }
                else if (!atoi(argv[i]))
                {
                    printf("\n****ERROR! Precision must be specified\n");
                    return SLIP_INCORRECT_INPUT;
                }
                option->prec = atoi(argv[i]);
                if (option->prec < 2)
                {
                    printf("\n\n****ERROR! Invalid precision. prec >= 2\n");
                    return SLIP_INCORRECT_INPUT;
                }
            }
        }
        else if ( strcmp(arg, "f") == 0 || strcmp(arg, "file") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! Matrix name must be entered\n");
                return SLIP_INCORRECT_INPUT;
            }
            *mat_name = argv[i];
            if (!argv[++i])
            {
                printf("\n****ERROR! Right hand side vector name must"
                    " be entered\n");
                return SLIP_INCORRECT_INPUT;
            }
            *rhs_name = argv[i];
        }
        else
        {
            printf("\n\n**ERROR! Unknown command line parameter: %s"
                    "\nIgnoring this parameter\n",arg);
            return SLIP_INCORRECT_INPUT;
        }
    }
    return SLIP_OK;
}


/* Purpose: This function shows the usage of the code.*/
void SLIP_show_usage() //display the usage of the code
{
    printf("\n\n\t\t****USAGE****"
    "\n\t./SLIP_LU followed by:"
    "\n\tc: Indicates soln will be checked\
     \n\tp (or piv) 0~5 : indicate type of pivoting"
    "\n\tcol or q: column order used: 0: none, 1: COLAMD, 2: AMD"
    "\n\tt or tol: tolerance parameter\
     \n\to2 or out2: output printed to screen"
    "\n\tf or file: filenames. must be of format MATRIX_NAME RHS_NAME"
    "\n\to or out: output will be printed to file. Must be followed by\
         1: rational, 2: double, 3 PREC: float of precision PREC"
    "\n****REFER TO README.txt FOR DETAILED DESCRIPTION OF INPUT PARAMETERS****\
     \n");
}


/* Purpose: This function reads in a matrix stored in a triplet format
 * This format used can be seen in any of the example mat files. 
 * 
 * This is only used for Demo purposes
 */
SLIP_info SLIP_tripread
(
    SLIP_sparse* A,     // Matrix to be populated
    FILE* file          // file to read from (must already be open)
)
{
    SLIP_info ok;
    if (A == NULL || file == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    int32_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    ok = fscanf(file, "%d %d %d\n", &m, &n, &nz);
    if (feof(file) || ok < 3)
    {
        return SLIP_INCORRECT_INPUT;
    }
    // Initialize i and j vectors 
    int32_t *i = (int32_t*) SLIP_malloc(nz * sizeof(int32_t));
    int32_t *j = (int32_t*) SLIP_malloc(nz * sizeof(int32_t));

    // Create an initialized input mpz vector
    mpz_t* x_mpz = SLIP_create_mpz_array(nz);
    
    if (!i || !j || !x_mpz)
    {
        SLIP_FREE(i);
        SLIP_FREE(j);
        SLIP_delete_mpz_array(&x_mpz, nz);
        return SLIP_OUT_OF_MEMORY;
    }

    int32_t decrement;
    ok = SLIP_gmp_fscanf(file, "%d %d %Zd\n", &i[0], &j[0], &x_mpz[0]);
    if (feof(file) || ok < 3)
    {
        SLIP_FREE(i);
        SLIP_FREE(j);
        SLIP_delete_mpz_array(&x_mpz, nz);
        return SLIP_INCORRECT_INPUT;
    }

    if (SLIP_MIN(i[0], j[0]) == 0)
    {
        decrement = 0;
    }
    else
    {
        decrement = 1;
        i[0]-=decrement;
        j[0]-=decrement;
    }

    // Read in the values from file
    for (int32_t p = 1; p < nz; p++)
    {
        ok = SLIP_gmp_fscanf(file, "%d %d %Zd\n", &i[p], &j[p], &x_mpz[p]); 
        if ((feof(file) && p != nz-1) || ok < 3)
        {
            SLIP_FREE(i);
            SLIP_FREE(j);
            SLIP_delete_mpz_array(&x_mpz, nz);
            return SLIP_INCORRECT_INPUT;
        }
        // Conversion from 1 based to 0 based if necessary
        i[p] -= decrement;
        j[p] -= decrement;
    }

    //------------------------------------------------------------------
    // At this point, we have read in i, j, and x arrays and have 
    // allocated memory for the A matrix. The i & j are stored as 
    // int32_t and x is stored as a mpz_array. We conclude by using the 
    // appropriate SLIP_build_* to construct our input matrix A
    //------------------------------------------------------------------
    ok = SLIP_build_sparse_trip_mpz(A, i, j, x_mpz, n, nz);
    
    // A now contains our input matrix. Free memory for i, j, and x
    
    SLIP_FREE(i);
    SLIP_FREE(j);
    SLIP_delete_mpz_array(&x_mpz, nz);
    return ok;
}


/* Purpose: This function reads in a double matrix stored in a triplet format
 * This format used can be seen in any of the example mat files. 
 * 
 * This is only used for Demo purposes
 */

SLIP_info SLIP_tripread_double
(
    SLIP_sparse* A,        // Matrix to be populated
    FILE* file,          // file to read from (must already be open)
    SLIP_options* option
)
{
    SLIP_info ok;
    if (A == NULL || file == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }
    // Read in triplet form first
    int32_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    ok = fscanf(file, "%d %d %d\n", &m, &n, &nz);
    if (feof(file) || ok < 3)
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    int32_t *i = (int32_t*) SLIP_malloc(nz* sizeof(int32_t));
    int32_t *j = (int32_t*) SLIP_malloc(nz* sizeof(int32_t));
    double *x_doub = (double*) SLIP_malloc(nz* sizeof(double));

    if (!i || !j || !x_doub )
    {
        SLIP_FREE(i);                 
        SLIP_FREE(j);                     
        SLIP_FREE(x_doub);                
        return SLIP_OUT_OF_MEMORY;
    }

    int32_t decrement;
    ok = fscanf(file, "%d %d %lf\n", &(i[0]), &(j[0]), &(x_doub[0]));
    if (feof(file) || ok < 3)
    {
        SLIP_FREE(i);                 
        SLIP_FREE(j);                     
        SLIP_FREE(x_doub);                
        return SLIP_INCORRECT_INPUT;
    }

    if (SLIP_MIN(i[0], j[0]) == 0)
    {
        decrement = 0;
    }
    else
    {
        decrement = 1;
        i[0]-=decrement;
        j[0]-=decrement;
    }

    // Read in the values from file
    for (int32_t k = 1; k < nz; k++)
    {
        ok = fscanf(file, "%d %d %lf\n", &(i[k]), &(j[k]), &(x_doub[k]));
        if ((feof(file) && k != nz-1) || ok < 3)
        {
            SLIP_FREE(i);                 
            SLIP_FREE(j);                     
            SLIP_FREE(x_doub);                
            return SLIP_INCORRECT_INPUT;
        }
        // Conversion from 1 based to 0 based
        i[k] -= decrement;
        j[k] -= decrement;
    }

    //------------------------------------------------------------------
    // At this point, we have read in i, j, and x arrays and have 
    // allocated memory for the A matrix. The i & j are stored as 
    // int32_t and x is stored as a double array. We conclude by using the
    // appropriate SLIP_build_* to construct our input matrix A
    //------------------------------------------------------------------
    
    ok = SLIP_build_sparse_trip_double(A, i, j, x_doub, n, nz, option);
    // Now, A contains the input matrix
    SLIP_FREE(i);                 
    SLIP_FREE(j);                     
    SLIP_FREE(x_doub);                
    return ok;
}


/* Purpose: Read a dense matrix. This is for demo purposes only */

SLIP_info SLIP_read_dense
(
    SLIP_dense *b,
    FILE* file          // file to read from (must already be open)
)
{
    
    //------------------------------------------------------------------
    // We read in a dense matrix and then utilize the appropriate
    // SLIP_build_dense_*. Here, we assume that the input is read in
    // as a dense mpz_t** matrix. The main component of this code is 
    // reading in said matrix.
    //------------------------------------------------------------------
    
    if (b == NULL || file == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t nrows, ncols;
    SLIP_info ok;

    // First, we obtain the dimension of the matrix
    ok = fscanf(file, "%d %d", &nrows, &ncols);
    if (feof(file) || ok < 2)
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    // Now, we create our dense mpz_t matrix
    mpz_t** b_orig = SLIP_create_mpz_mat(nrows, ncols);
    if (b_orig == NULL)
    {
        return SLIP_OUT_OF_MEMORY;
    }
    
    // We now populate the matrix b.
    for (int32_t i = 0; i < nrows; i++)
    {
        for (int32_t j = 0; j < ncols; j++)
        {
            ok = SLIP_gmp_fscanf(file, "%Zd", &(b_orig[i][j])); 
            if (ok < 0) // return from this function can be a nonzero
            {
                        printf("\n\nhere at i = %d and j = %d", i, j);

                return SLIP_INCORRECT_INPUT;
            }
        }
    }
    
    //------------------------------------------------------------------
    // At this point, b_orig contains our original matrix stored in
    // mpz_t** format. We now utilize the appropriate SLIP_build function
    // to form our internal SLIP_dense structure.
    //------------------------------------------------------------------
    
    ok = SLIP_build_dense_mpz(b, b_orig, nrows, ncols);
    
    return SLIP_OK;
}

#define SLIP_PRINT(...)      \
{                            \
    if (out_file == NULL )   \
    {                        \
        printf(__VA_ARGS__); \
    }                        \
    else                     \
    {                        \
        fprintf(out_file, __VA_ARGS__); \
    }                        \
}

//------------------------------------------------------------------------------
// SLIP_print_stats_mpq: prints the solution vector(s) as a set of mpq_t entries
//------------------------------------------------------------------------------

SLIP_info SLIP_print_stats_mpq
(
    FILE *out_file,         // file to print to
    mpq_t **x_mpq,          // solution vector in mpq, pass NULL if unused
    int32_t n,              // dimension of A
    int32_t numRHS,         // number of RHS vectors
    SLIP_info check,        // whether the solution is correct or not
    SLIP_options *option    // option struct telling how much info to print
)
{
    SLIP_info ok = SLIP_OK;
    if (option == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

    // Info about output file
    if (option->print_level >= 2 && out_file != NULL)
    {
        if (x_mpq == NULL)
        {
            return SLIP_INCORRECT_INPUT;
        }
        fprintf(out_file,
            "\nSolution output in full precision rational arithmetic\n");
        for (int32_t i = 0; i < n; i++)
        {
            for (int32_t j = 0; j < numRHS; j++)
            {
                ok = SLIP_gmp_fprintf(out_file, "%Qd ", x_mpq[i][j]); 
                if (ok < 0)
                {
                    return ok;
                }
            }
            fprintf(out_file, "\n");
        }
    }
    return SLIP_OK;
}

//------------------------------------------------------------------------------
// SLIP_print_stats_double:  prints the solution vector(s) as a set of double entries
//------------------------------------------------------------------------------

SLIP_info SLIP_print_stats_double
(
    FILE *out_file,         // file to print to
    double **x_doub,        // solution vector in double, pass NULL if unused
    int32_t n,              // dimension of A
    int32_t numRHS,         // number of RHS vectors
    SLIP_info check,        // whether the solution is correct or not
    SLIP_options *option    // option struct telling how much info to print
)
{
    if (option == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

    // Info about output file
    if (option->print_level >= 2 && out_file != NULL)
    {
        if (x_doub == NULL)
        {
            return SLIP_INCORRECT_INPUT;
        }
        fprintf(out_file, "\nSolution output in double precision\n");
        for (int32_t i = 0; i < n; i++)
        {
            for (int32_t j = 0; j < numRHS; j++)
            {
                // Output the solution in double precision
                fprintf(out_file, "%lf ", x_doub[i][j]);
            }
            fprintf(out_file,"\n");
        }
    }
    return SLIP_OK;
}

//------------------------------------------------------------------------------
// SLIP_print_stats_mpfr:prints the solution vector(s) as a set of mpfr_t entries
//------------------------------------------------------------------------------

SLIP_info SLIP_print_stats_mpfr
(
    FILE *out_file,         // file to print to
    mpfr_t **x_mpfr,        // solution vector in mpfr, pass NULL if unused
    int32_t n,              // dimension of A
    int32_t numRHS,         // number of RHS vectors
    SLIP_info check,        // whether the solution is correct or not
    SLIP_options *option    // option struct telling how much info to print
)
{
    SLIP_info ok = SLIP_OK;
    if (option == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

   // Info about output file
    if (option->print_level >= 2 && out_file != NULL)
    {
        if (x_mpfr == NULL)
        {
            return SLIP_INCORRECT_INPUT;
        }
        fprintf(out_file, "\nSolution output in fixed precision of size:"
            " %ld bits\n", option->prec);
        for (int32_t i = 0; i < n; i++)
        {
            for (int32_t j = 0; j < numRHS; j++)
            {
                ok = SLIP_mpfr_fprintf(out_file, "%.*Rf",
                    option->prec, x_mpfr[i][j]); 
                if (ok < 0)
                {
                    return ok;
                }
            }
            fprintf(out_file, "\n");
        }
    }
    return SLIP_OK;
}
