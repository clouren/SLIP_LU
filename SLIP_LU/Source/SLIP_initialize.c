//------------------------------------------------------------------------------
// SLIP_LU/SLIP_initialize: initialize SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// SLIP_initialize initializes the working evironment for SLIP_LU.

#include "slip_internal.h"

void SLIP_initialize ( void )
{
    mp_set_memory_functions (slip_gmp_allocate, slip_gmp_reallocate,
        slip_gmp_free) ;
}

