//------------------------------------------------------------------------------
// SLIP_LU/Demo/demos.h: #include file the demo programs
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU.h"
#include <time.h>
#include <stdint.h>
#include <inttypes.h>

#define OK(method)                      \
{                                       \
    ok = method ;                       \
    if (ok != SLIP_OK)                  \
    {                                   \
        printf ("Error: line %d file %s\n", __LINE__, __FILE__) ; \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}

#define SLIP_MIN(a,b) (((a) < (b)) ? (a) : (b))

/* Purpose: This processes the command line for user specified options */
SLIP_info SLIP_process_command_line //processes the command line
(
    int argc,               // number of command line arguments
    char* argv[],           // set of command line arguments
    SLIP_options* option,   // struct containing the command options
    char** mat_name,        // Name of the matrix to be read in
    char** rhs_name,        // Name of the RHS vector to be read in
    SLIP_type *rat          // data type of output solution:
                            // 1:SLIP_MPZ (default), 2:SLIP_FP64, 3:SLIP_MPFR
);

/* Purpose: This function prints out the user specified/default options*/
void SLIP_print_options     // display specified/default options to user
(
    SLIP_options* option // struct containing all of the options
);

/* Purpose: This function shows the usage of the code.*/
void SLIP_show_usage(void);

/* Purpose: This function reads in a matrix stored in a triplet format.
 * This format used can be seen in any of the example mat files.
 */
SLIP_info SLIP_tripread
(
    SLIP_sparse **A_handle,     // Matrix to be constructed
    FILE* file                  // file to read from (must already be open)
) ;

/* Purpose: This function reads in a double matrix stored in a triplet format.
 * This format used can be seen in any of the example mat files.
 */
SLIP_info SLIP_tripread_double
(
    SLIP_sparse **A_handle,     // Matrix to be constructed
    FILE* file,                 // file to read from (must already be open)
    SLIP_options* option
) ;

/* Purpose: SLIP_read_dense: read a dense matrix. */
SLIP_info SLIP_read_dense
(
    SLIP_dense **b_handle,      // Matrix to be constructed
    FILE* file                  // file to read from (must already be open)
) ;

/* Purpose: prints the solution vector(s) as a set of mpq_t entries */
SLIP_info SLIP_print_stats_mpq
(
    mpq_t **x_mpq,          // solution vector in mpq, pass NULL if unused
    int64_t n,              // dimension of A
    int64_t numRHS,         // number of RHS vectors
    SLIP_info check,        // whether the solution is correct or not
    SLIP_options *option    // option struct telling how much info to print
);

/* Purpose: prints the solution vector(s) as a set of double entries */
SLIP_info SLIP_print_stats_double
(
    double **x_doub,        // solution vector in double, pass NULL if unused
    int64_t n,              // dimension of A
    int64_t numRHS,         // number of RHS vectors
    SLIP_info check,        // whether the solution is correct or not
    SLIP_options *option    // option struct telling how much info to print
);

/* Purpose: prints the solution vector(s) as a set of mpfr_t entries */
SLIP_info SLIP_print_stats_mpfr
(
    mpfr_t **x_mpfr,        // solution vector in mpfr, pass NULL if unused
    int64_t n,              // dimension of A
    int64_t numRHS,         // number of RHS vectors
    SLIP_info check,        // whether the solution is correct or not
    SLIP_options *option    // option struct telling how much info to print
);

