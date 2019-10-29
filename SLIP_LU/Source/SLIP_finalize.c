# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function initialize the working evironment for SLIP LU library.
 */
void SLIP_finalize
(
    void
)
{
    slip_mpfr_free_cache();       // Free mpfr internal cache
    slip_gmp_finalize();          // Reset GMP memory variables 
}
