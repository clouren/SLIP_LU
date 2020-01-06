//------------------------------------------------------------------------------
// SLIP_LU/SLIP_initialize_expert: intialize SLIP_LU memory functions for GMP
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

# include "SLIP_LU_internal.h"

/* 
 * Purpose: This function initializes the working environment for SLIP_LU with
 * custom memory functions that are used for GMP. If the user passes in their
 * own malloc, realloc, or free function(s) we use those internally to process
 * memory. If a NULL pointer is passed in for any function, then default
 * functions are used.
 *
 * The three functions are similar to ANSI C malloc, realloc, and free
 * functions, but the calling syntax is not the same.  Below are the defintions
 * that must be followed, per the GMP specification:
 *
 *  void *MyMalloc (size_t size) ;  // same as the ANSI C malloc
 *  void *MyRealloc (void *p, size_t oldsize, size_t newsize) ; // differs
 *  void MyFree (void *p, size_t size) ; // differs
 *
 * MyMalloc has the identical to the ANSI C malloc  MyRealloc adds a parameter,
 * oldsize, which is the prior size of the block of memory to be reallocated.
 * MyFree takes a second argument, which is the size of the block that is being
 * freed.
 *
 * The default memory management functions used inside GMP are:
 *
 *      MyMalloc    slip_gmp_allocate
 *      MyRealloc   slip_gmp_reallocate
 *      MyFree      slip_gmp_free
 *
 * The SLIP_gmp_* memory management functions are unique to SLIP_LU.  They
 * provide an elegant workaround for how GMP manages its memory.  By default,
 * if GMP attempts to allocate memory, but it fails, then it simply terminates
 * the user application.  This behavoir is not suitable for many applications
 * (MATLAB in particular).  Fortunately, GMP allows the user application
 * (SLIP_LU in this case) to pass in alternative memroy manager functions, via
 * mp_set_memory_functions.  The SLIP_gmp_* functions do not return to GMP if
 * the allocation fails, but instead use the longjmp feature of ANSI C to
 * implement a try/catch mechanism.  The memory failure can then be safely
 * handled by SLIP_LU, without memory leaks and without terminating the user
 * application.
 *
 * When SLIP_LU is used via MATLAB, the following functions are used instead:
 *
 *      MyMalloc    mxMalloc
 *      MyRealloc   SLIP_gmp_mex_realloc (a wrapper for mxRealloc)
 *      MyFree      SLIP_gmp_mex_free (a wrapper for mxFree)
 *
 * Note that these functions are not used by SLIP_LU itself, but only inside
 * GMP.  The functions used by SLIP_LU itself are SLIP_malloc, SLIP_calloc,
 * SLIP_realloc, and SLIP_free, which are wrappers for the ANSI C malloc,
 * calloc, realloc, and free, or (if used inside MATLAB), for the MATLAB
 * mxMalloc, mxCalloc, mxRealloc, and mxFree functions.
 */

void SLIP_initialize_expert
(
    void* (*MyMalloc) (size_t),                     // User defined malloc function
    void* (*MyRealloc) (void *, size_t, size_t),    // User defined realloc function
    void (*MyFree) (void*, size_t)                  // User defined free function
)
{

    //--------------------------------------------------------------------------
    // Set defaults for any NULL pointers passed in
    //--------------------------------------------------------------------------

    if (MyMalloc == NULL)
    {
        MyMalloc = slip_gmp_allocate ;
    }

    if (MyRealloc == NULL)
    {
        MyRealloc = slip_gmp_reallocate ;
    }

    if (MyFree == NULL)
    {
        MyFree = slip_gmp_free ;
    }

    //--------------------------------------------------------------------------
    // Set GMP memory functions 
    //--------------------------------------------------------------------------

    mp_set_memory_functions ((*MyMalloc), (*MyRealloc), (*MyFree)) ;
}

