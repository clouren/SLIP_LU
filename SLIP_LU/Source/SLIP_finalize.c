//------------------------------------------------------------------------------
// SLIP_LU/SLIP_finalize: finalize SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// SLIP_finalize frees the working environment for SLIP LU library.

#include "slip_internal.h"

void SLIP_finalize
(
    void
)
{
    SLIP_mpfr_free_cache ( ) ;    // Free mpfr internal cache
    slip_gmp_finalize ( ) ;       // Reset GMP memory variables
}

