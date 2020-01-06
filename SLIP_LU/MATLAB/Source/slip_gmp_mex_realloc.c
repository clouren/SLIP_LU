//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/SLIP_gmp_mex_realloc: A gmp realloc function for matlab
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"

/* Purpose: A GMP reallocation function 
 * This allows GMP to use MATLAB's default realloc function 
 */

// A GMP realloc function
void* SLIP_gmp_mex_realloc 
(
    void* x,    // void* to be reallocated 
    size_t a,   // Previous size
    size_t b    // New size
)
{
    return mxRealloc(x,b);
}
