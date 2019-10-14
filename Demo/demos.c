# include "demos.h"

/* Purpose: This function prints out the user specified/default options*/
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
        order = "the UMFPACK";
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
        else if ( strcmp(arg,"c") == 0 || strcmp(arg,"check") == 0)
        {
     	    option->check = true;
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
