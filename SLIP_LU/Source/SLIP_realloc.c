//------------------------------------------------------------------------------
// SLIP_LU/SLIP_realloc: wrapper for realloc
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

// A wrapper for realloc

// If p is non-NULL on input, it points to a previously allocated object of
// size old_size * size_of_item.  The object is reallocated to be of size
// new_size * size_of_item.  If p is NULL on input, then a new object of that
// size is allocated.  On success, a pointer to the new object is returned, and
// ok is returned as true.  If the allocation fails, ok is set to false and a
// pointer to the old (unmodified) object is returned.

void* SLIP_realloc 
(
    void* p,            // Pointer to be realloced
    size_t old_size,    // Old size of this pointer
    size_t new_size     // New size of this pointer
)
{
    #ifdef MATLAB_MEX_FILE
        return mxRealloc(p, new_size);
    #else
        // Ensure at least one byte is allocated
        old_size = SLIP_MAX(1, old_size);
        new_size = SLIP_MAX(1, new_size);
        
        if (p == NULL)
        {
            return SLIP_malloc(new_size);
        }
        else if (new_size == old_size)
        {
            return p;
        }
        else
        {
            void* pnew;
            pnew = (void*) SLIP_MEMORY_REALLOC (p, new_size) ;
            if (pnew == NULL)
            {
                if (new_size < old_size)
                {
                    // The attempt to reduce the size of the block failed,
                    // but the old block is unchanged. Pretend to succeed
                    return p;
                }
                else
                {
                    // Out of memory
                    SLIP_FREE(p);
                    return NULL;
                }
            }
            else
            {
                p = pnew;
                return p;
            }
        }
    #endif
}

