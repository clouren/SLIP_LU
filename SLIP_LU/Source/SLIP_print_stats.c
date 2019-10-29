# include "SLIP_LU_internal.h"

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

// print the correctness of the solution if option->check enabled
static inline SLIP_info print_check
(
    FILE *out_file,         // file to print to
    SLIP_info check,        // whether the solution is correct or not
    SLIP_options *option    // option struct telling how much info to print
)
{
    // Was the solution checked?
    if (option->check)
    {
        if (check == SLIP_OK)
        {
            if (option->print_level > 0)
            {
                SLIP_PRINT("\nSolution is correct!\n");
            }
        }

        // Incorrect solution. This should not happen
        else
        {
            if (option->print_level > 0)
            {
                SLIP_PRINT("\n****ERROR**** Solution incorrect!"
                       "\n\nHave you modified the source code?"
                       "\nReinstall and if this issue persists email me at:"
                       " clouren@tamu.edu\n\n");
            }
            return SLIP_INCORRECT;
        }
    }
    return SLIP_OK;
}

/* ========================================================================== */
/* This function prints statistics about the SLIP LU factorization            */
/* ========================================================================== */
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
    // Is the solution correct?
    if (print_check (out_file, check, option) != SLIP_OK)
    {
        return SLIP_INCORRECT;
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
                ok = slip_gmp_fprintf(out_file, "%Qd ", x_mpq[i][j]);
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
    // Is the solution correct?
    if (print_check (out_file, check, option) != SLIP_OK)
    {
        return SLIP_INCORRECT;
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
    // Is the solution correct?
    if (print_check (out_file, check, option) != SLIP_OK)
    {
        return SLIP_INCORRECT;
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
                ok = slip_mpfr_fprintf(out_file, "%.*Rf",
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
#undef SLIP_PRINT
