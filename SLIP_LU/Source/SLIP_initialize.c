//------------------------------------------------------------------------------
// SLIP_LU/SLIP_initialize: initialize SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function initializes the working evironment for SLIP_LU.
 */

#include "slip_internal.h"

void SLIP_initialize ( void )
{

    //--------------------------------------------------------------------------
    // Set GMP memory functions as default SLIP gmp functions
    //--------------------------------------------------------------------------

    SLIP_initialize_expert (NULL, NULL, NULL) ;
}

