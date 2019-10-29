//------------------------------------------------------------------------------
// SLIP_LU/SLIP_initialize: intialize SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function initializes the working evironment for SLIP_LU
 */
void SLIP_initialize( void )
{
    //--------------------------------------------------------------------------
    // Set GMP memory functions 
    //--------------------------------------------------------------------------

    mp_set_memory_functions(slip_gmp_allocate, slip_gmp_reallocate,
            slip_gmp_free);
}

