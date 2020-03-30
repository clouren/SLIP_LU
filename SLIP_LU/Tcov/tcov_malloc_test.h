//------------------------------------------------------------------------------
// SLIP_LU/Tcov/tcov_malloc_test.h
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#ifndef SLIP_TCOV_MALLOC_TEST_H
#define SLIP_TCOV_MALLOC_TEST_H

#include "SLIP_LU_internal.h"

#ifdef SLIP_MEMORY_REALLOC
#undef SLIP_MEMORY_REALLOC
#endif
void* slip_realloc_wrapper
(
    void* p,           // Pointer to be realloced
    size_t new_size    // Size to alloc
);
/* to be used in SLIP_gmp.c */
#define SLIP_MEMORY_REALLOC slip_realloc_wrapper

extern int64_t malloc_count;
#define GOTCHA \
    printf ("%s, line %d, slip_gmp_ntrials = %ld, malloc_count = %ld\n", \
    __FILE__, __LINE__, slip_gmp_ntrials, malloc_count);

#define SLIP_PRINT_INFO(info)                                               \
{                                                                           \
    printf ("file %s line %d: ", __FILE__, __LINE__) ;                      \
    switch(info)                                                            \
    {                                                                       \
        case SLIP_OK:              printf("SLIP_OK\n");            break;   \
        case SLIP_OUT_OF_MEMORY:   printf("OUT OF MEMORY\n");      break;   \
        case SLIP_SINGULAR:        printf("Matrix is SINGULAR\n"); break;   \
        case SLIP_INCORRECT_INPUT: printf("INCORRECT INPUT\n");    break;   \
        case SLIP_INCORRECT:       printf("SLIP_INCORRECT\n");     break;   \
        default:                   printf("unknown!\n");                    \
    }                                                                       \
}

#ifdef SLIP_CHECK
#undef SLIP_CHECK
#endif

#define SLIP_CHECK(method)          \
{                                   \
    info = (method) ;               \
    if (info != SLIP_OK)            \
    {                               \
        SLIP_PRINT_INFO (info)      \
        SLIP_FREE_ALL ;             \
        return (info) ;             \
    }                               \
}

SLIP_info slip_gmp_realloc_test
(
    void **p_new,
    void * p_old,
    size_t old_size,
    size_t new_size
) ;

#endif
