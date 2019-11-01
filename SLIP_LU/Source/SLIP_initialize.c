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
    // Set GMP memory functions as default SLIP gmp functions
    //--------------------------------------------------------------------------

    SLIP_initialize_expert(NULL, NULL, NULL);
}

