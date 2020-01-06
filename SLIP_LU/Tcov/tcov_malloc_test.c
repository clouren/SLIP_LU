#include "tcov_malloc_test.h"
//#include "SLIP_LU_internal.h"
int malloc_count = INT_MAX ;

// Note that only the ANSI C memory manager is used here
// (malloc, calloc, realloc, free)

/* wrapper for malloc */
void* SLIP_malloc
(
    size_t size        // Size to alloc
)
{
    if (--malloc_count < 0) {printf("malloc pretend to fail\n");return (NULL) ;} /* pretend to fail */
    // ensure at least one byte is malloc'd
    if (size <= 0) {size = 1 ;}
    return (malloc(size));
}

/* wrapper for calloc */
void* SLIP_calloc
(
    size_t n,          // Size of array
    size_t size        // Size to alloc
)
{
    if (--malloc_count < 0) {printf("calloc pretend to fail\n");return (NULL) ;} /* pretend to fail */
    // ensure at least one byte is calloc'd
    if (n <= 0) {n = 1 ;}
    if (size <= 0) {size = 1 ;}
    return (calloc(n,size));
}

void* slip_realloc_wrapper
(
    void* p,           // Pointer to be realloced
    size_t new_size    // Size to alloc
)
{
    if (--malloc_count < 0)
    {
        /* pretend to fail */
        printf("realloc pretend to fail\n");
        return (NULL);
    }
    else
    {
        return (realloc(p, new_size));
    }
}

/* wrapper for realloc */
void* SLIP_realloc 
(
    void* p,            // Pointer to be realloced
    size_t old_size,    // Old size of this pointer
    size_t new_size     // New size of this pointer
)
{
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
        void* pnew = (void*) SLIP_MEMORY_REALLOC (p, new_size);
        if (pnew == NULL)
        {
            if (new_size <= old_size)
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
}

jmp_buf SLIP_gmp_environment ;  // for setjmp and longjmp

int SLIP_gmp_realloc_test
(
    void **p_new,
    void * p_old,
    size_t old_size,
    size_t new_size
)
{
    int SLIP_gmp_status = setjmp (SLIP_gmp_environment);
    if (SLIP_gmp_status != 0)
    {
        return SLIP_OUT_OF_MEMORY;
    }
    *p_new = SLIP_gmp_reallocate(p_old, old_size, new_size);
    return SLIP_OK;
}
