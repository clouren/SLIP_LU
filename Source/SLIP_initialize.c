# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function initialize the working evironment for SLIP LU library.
 */
void SLIP_initialize( void )
{
    //--------------------------------------------------------------------------
    // Set GMP memory functions 
    //--------------------------------------------------------------------------

    mp_set_memory_functions(slip_gmp_allocate, slip_gmp_reallocate,
            slip_gmp_free);
}
