#include "SLIP_LU_internal.h"


#define OK(method)                      \
{                                       \
    ok = method ;                       \
    if (ok != SLIP_OK)                  \
    {                                   \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
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

);
 
/* Purpose: This function prints out the user specified/default options*/ 
void SLIP_print_options     // display specified/default options to user 
( 
    SLIP_options* option // struct containing all of the options 
); 

/* Purpose: This function shows the usage of the code.*/
void SLIP_show_usage(void);

