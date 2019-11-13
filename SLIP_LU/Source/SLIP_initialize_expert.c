//------------------------------------------------------------------------------
// SLIP_LU/SLIP_initialize: intialize SLIP_LU with user defined memory functions
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function initializes the working evironment for SLIP_LU
 */
void SLIP_initialize_expert
(
    void* (*MyMalloc) (size_t),                     // User defined malloc function
    void* (*MyRealloc) (void *, size_t, size_t),    // User defined realloc function
    void (*MyFree) (void*, size_t)                  // User defined free function
)
{
    //--------------------------------------------------------------------------
    // Set GMP memory functions 
    //--------------------------------------------------------------------------

    //MyMalloc = slip_gmp_allocate;
    //MyRealloc = slip_gmp_reallocate;
    //MyFree = slip_gmp_free;
    
    //TODO make this temporary/fix?
    if (MyMalloc == 0)
        MyMalloc = slip_gmp_allocate;
    if (MyRealloc == 0)
        MyRealloc = slip_gmp_reallocate;
    if (MyFree == 0)
        MyFree = slip_gmp_free;

    mp_set_memory_functions((*MyMalloc), (*MyRealloc),
            (*MyFree));
}

