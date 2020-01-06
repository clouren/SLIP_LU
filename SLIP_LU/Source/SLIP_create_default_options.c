//------------------------------------------------------------------------------
// SLIP_LU/SLIP_create_default_options: set defaults
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* Purpose: Create and return SLIP_options pointer with default parameters
 * upon successful allocation, which are defined in SLIP_LU_internal.h
 */

SLIP_options* SLIP_create_default_options ( void )
{
    SLIP_options* option = SLIP_malloc(sizeof(SLIP_options));
    if (!option) {return NULL;}
    option->pivot    = SLIP_DEFAULT_PIVOT;
    option->order    = SLIP_DEFAULT_ORDER;
    option->print_level= SLIP_DEFAULT_PRINT_LEVEL;
    option->prec     = SLIP_DEFAULT_PRECISION;
    option->tol      = SLIP_DEFAULT_TOL;
    option->SLIP_MPFR_ROUND = SLIP_DEFAULT_MPFR_ROUND;
    return option;
}

