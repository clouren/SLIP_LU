//------------------------------------------------------------------------------
// SLIP_LU/SLIP_gmp: interface to the gmp library
//------------------------------------------------------------------------------


// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// This file (SLIP_gmp.c) provides a wrapper for all functions in the GMP
// library used by SLIP_LU.  The wrappers enable memory failures to be caught
// and handled properly.  GMP, by default, aborts the user's application if any
// internal malloc fails.  This is not acceptable in a robust end-user
// application.  Fortunately, GMP allows the user package (SLIP_LU in this
// case) to pass in function pointers for malloc, calloc, realloc, and free.
// These functions are defined below.  If they fail, they do not return to GMP.
// Instead, they use the ANSI C longjmp feature to trap the error, and return
// the error condition to the caller.

# include "SLIP_LU_internal.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

jmp_buf slip_gmp_environment ;  // for setjmp and longjmp
int64_t slip_gmp_nmalloc = 0 ;  // number of malloc'd objects in SLIP_gmp_list
int64_t slip_gmp_nlist = 0 ;    // size of the SLIP_gmp_list
void **slip_gmp_list = NULL ;   // list of malloc'd objects

int64_t slip_gmp_ntrials = -1 ; // number of malloc's allowed (for
                                // testing only): -1 means unlimited.

mpz_t  *slip_gmpz_archive  = NULL ;    // current mpz object
mpq_t  *slip_gmpq_archive  = NULL ;    // current mpq object
mpfr_t *slip_gmpfr_archive = NULL ;    // current mpfr object

//------------------------------------------------------------------------------
// slip_gmp_init: initialize gmp
//------------------------------------------------------------------------------

/* Purpose: Create the list of malloc'd objects. This should be called before
 * calling any GMP function. It is also called by SLIP_gmp_allocate when
 * SLIP_gmp_list is NULL
 */

bool slip_gmp_init ( )
{
    slip_gmp_nmalloc = 0 ;
    slip_gmp_nlist = SLIP_GMP_LIST_INIT ;
    slip_gmp_list = (void**) SLIP_malloc (slip_gmp_nlist* sizeof (void *)) ;
    return (slip_gmp_list != NULL) ;
}

//------------------------------------------------------------------------------
// SLIP_gmp_finalize: finalize gmp
//------------------------------------------------------------------------------

/* Purpose: Free the list. Must be called when all use of GMP is done */
void slip_gmp_finalize ( )
{
    slip_gmpz_archive = NULL ;
    slip_gmpq_archive = NULL ;
    slip_gmpfr_archive = NULL ;
    slip_gmp_nmalloc = 0 ;
    slip_gmp_nlist = 0 ;
    SLIP_FREE (slip_gmp_list) ;
}

//------------------------------------------------------------------------------
// SLIP_gmp_allocate: malloc space for gmp
//------------------------------------------------------------------------------

/* Purpose: malloc space for gmp. A NULL pointer is never returned to the GMP
 * library. If the allocation fails, all memory allocated since the start of
 * the SLIP_gmp_wrapper is freed and an error is thrown to the GMP wrapper via
 * longjmp
 */

void *slip_gmp_allocate
(
    size_t size // Amount of memory to be allocated
)
{

    #ifdef SLIP_GMP_MEMORY_DEBUG
    printf ("slip_gmp_malloc (%g): ", (double) size) ;
    #endif

    //--------------------------------------------------------------------------
    // for testing only:
    //--------------------------------------------------------------------------

    if (slip_gmp_ntrials == 0)
    {
        // pretend to fail
        #ifdef SLIP_GMP_MEMORY_DEBUG
        printf ("slip_gmp_allocate pretends to fail\n") ;
        #endif
        longjmp (slip_gmp_environment, 1) ;
    }
    else if (slip_gmp_ntrials > 0)
    {
        // one more malloc has been used up
        slip_gmp_ntrials-- ;
    }

    //--------------------------------------------------------------------------
    // ensure the SLIP_gmp_list is large enough
    //--------------------------------------------------------------------------

    if (slip_gmp_list == NULL)
    {
        // create the initial SLIP_gmp_list
        if (!slip_gmp_init ( ))
        {
            // failure to create the SLIP_gmp_list
            longjmp (slip_gmp_environment, 2) ;
        }
    }
    else if (slip_gmp_nmalloc == slip_gmp_nlist)
    {
        // double the size of the SLIP_gmp_list
        size_t newsize = 2 * slip_gmp_nlist * sizeof (void *) ;
        // cannot use SLIP_realloc here, since it frees its input on failure.
        void **newlist = (void**) SLIP_MEMORY_REALLOC (slip_gmp_list, newsize) ;

        if (newlist == NULL)
        {
            // failure to double the size of the SLIP_gmp_list.
            // The existing SLIP_gmp_list is still valid, with the old size,
            // (SLIP_gmp_nlist).  This is required so that the error handler
            // can traverse the SLIP_gmp_list to free all objects there.
            longjmp (slip_gmp_environment, 3) ;
        }

        // success; the old SLIP_gmp_list has been freed, and replaced with
        // the larger newlist.
        slip_gmp_list = newlist ;
        slip_gmp_nlist *= 2 ;
    }

    //--------------------------------------------------------------------------
    // malloc the block
    //--------------------------------------------------------------------------

    void *p = SLIP_malloc (size) ;

    if (p == NULL)
    {
        // failure to allocate the new block
        longjmp (slip_gmp_environment, 4) ;
    }

    //--------------------------------------------------------------------------
    // save p in the SLIP_gmp_list and return result to GMP
    //--------------------------------------------------------------------------

    slip_gmp_list [slip_gmp_nmalloc++] = p ;

    #ifdef SLIP_GMP_MEMORY_DEBUG
    printf (" %p\n", p) ;
    slip_gmp_dump ( ) ;
    #endif

    // return p to SLIP_gmp_function (NEVER return a NULL pointer to GMP!)
    assert (p != NULL) ;
    return (p) ;
}

//------------------------------------------------------------------------------
// slip_gmp_free: free space for gmp
//------------------------------------------------------------------------------

/* Purpose: Free space for GMP */
void slip_gmp_free
(
    void *p,        // Block to be freed
    size_t size     // Size of p
)
{
    #ifdef SLIP_GMP_MEMORY_DEBUG
    printf ("\n=================== free %p\n", p) ;
    slip_gmp_dump ( ) ;
    #endif

    if (p != NULL && slip_gmp_list != NULL)
    {
        // remove p from the SLIP_gmp_list
        for (int64_t i = 0 ; i < slip_gmp_nmalloc ; i++)
        {
            if (slip_gmp_list [i] == p)
            {
                #ifdef SLIP_GMP_MEMORY_DEBUG
                printf ("    found at i = %d\n", i) ;
                #endif
                slip_gmp_list [i] = slip_gmp_list [--slip_gmp_nmalloc] ;
                break ;
            }
        }
    }

    #ifdef SLIP_GMP_MEMORY_DEBUG
    slip_gmp_dump ( ) ;
    #endif

    // free p, even if it is not found in the SLIP_gmp_list.  p is only in the
    // SLIP_gmp_list if it was allocated inside the current GMP function.
    // If the block was allocated by one GMP function and freed by another,
    // it is not in the list.
    SLIP_SAFE_FREE (p) ;
}

//------------------------------------------------------------------------------
// slip_gmp_reallocate:  wrapper for realloc
//------------------------------------------------------------------------------

/* Purpose: Wrapper for GMP to call reallocation */
void *slip_gmp_reallocate
(
    void *p_old,        // Pointer to be realloc'd
    size_t old_size,    // Old size of p
    size_t new_size     // New size of p
)
{
    #ifdef SLIP_GMP_MEMORY_DEBUG
    printf ("slip_gmp_realloc (%p, %g, %g)\n", p_old,
        (double) old_size, (double) new_size) ;
    #endif

    if (p_old == NULL)
    {
        // realloc (NULL, size) is the same as malloc (size)
        return (slip_gmp_allocate (new_size)) ;
    }
    else if (new_size == 0)
    {
        // realloc (p, 0) is the same as free (p), and returns NULL
        slip_gmp_free (p_old, old_size) ;
        return (NULL) ;
    }
    else
    {
        // change the size of the block
        void *p_new = slip_gmp_allocate (new_size) ;
        // Note that p_new will never be NULL here, since SLIP_gmp_allocate
        // does not return if it fails.
        memcpy (p_new, p_old, SLIP_MIN (old_size, new_size)) ;
        slip_gmp_free (p_old, old_size) ;
        return (p_new) ;
    }
}

//------------------------------------------------------------------------------
// slip_gmp_dump: debug function
//------------------------------------------------------------------------------

/* Purpose: Dump the list of malloc'd objects */
#ifdef SLIP_GMP_MEMORY_DEBUG
void slip_gmp_dump ( )
{
    // dump the SLIP_gmp_list
    printf ("nmalloc = %g, SLIP_gmp_nlist = %g\n",
        (double) slip_gmp_nmalloc, (double) slip_gmp_nlist) ;
    if (slip_gmp_list != NULL)
    {
        for (int64_t i = 0 ; i < slip_gmp_nmalloc ; i++)
        {
            printf ("    slip_gmp_list [%d] = %p\n", i, slip_gmp_list [i]) ;
        }
    }
}
#endif

//------------------------------------------------------------------------------
// slip_gmp_failure: catch an error
//------------------------------------------------------------------------------

/* Purpose: Catch an error from longjmp */
void slip_gmp_failure
(
    int32_t status      // Status returned from longjmp
)
{
    #ifdef SLIP_GMP_MEMORY_DEBUG
    printf ("failure from longjmp: status: %d\n", status) ;
    #endif

    // first free all caches
    mpfr_free_cache ( ) ;

    // Free the list
    if (slip_gmp_list != NULL)
    {
        for (int64_t i = 0 ; i < slip_gmp_nmalloc ; i++)
        {
            SLIP_SAFE_FREE (slip_gmp_list [i]) ;
        }
    }
    slip_gmp_finalize ( ) ;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------------------Print and Scan functions--------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/* Safely print to the stream fp. Return positive value (the number of
 * characters written) upon success, otherwise return negative value (error
 * code) */
int32_t SLIP_gmp_fprintf(FILE *fp, const char *format, ... )
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call gmp_vfprintf
    int32_t num_of_char;
    va_list args;
    va_start(args, format);
    num_of_char = gmp_vfprintf(fp, format, args);
    va_end (args);

    // Finish the wrapper and return num_of_char if successful
    SLIP_GMP_WRAPPER_FINISH;
    // gmp_vfprintf returns -1 if an error occurred.
    return ((num_of_char < 0) ?
        ((int32_t) SLIP_INCORRECT_INPUT) : num_of_char) ;
}

/* Safely print to the standard output stdout. Return positive value (the number
 * of characters written) upon success, otherwise return negative value (error
 * code) */
int32_t SLIP_gmp_printf(const char *format, ... )
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call gmp_vprintf
    int32_t num_of_char;
    va_list args;
    va_start(args, format);
    num_of_char = gmp_vprintf(format, args);
    va_end (args);

    // Finish the wrapper and return num_of_char if successful
    SLIP_GMP_WRAPPER_FINISH;
    // gmp_vprintf returns -1 if an error occurred.
    return ((num_of_char < 0) ?
        ((int32_t) SLIP_INCORRECT_INPUT) : num_of_char) ;
}

/* Safely scan the stream fp. Return positive value (the number of fields
 * successfully parsed and stored), otherwise return negative value (error
 * code) */
int32_t SLIP_gmp_fscanf(FILE *fp, const char *format, ... )
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call gmp_vfscanf
    int32_t num_of_field;
    va_list args;
    va_start(args, format);
    num_of_field = gmp_vfscanf(fp, format, args);
    va_end (args);

    // Finish the wrapper and return num_of_field if successful
    SLIP_GMP_WRAPPER_FINISH;
    // If end of input (or a file error) is reached before a character
    // for a field or a literal, and if no previous non-suppressed fields have
    // matched, then the return value is EOF instead of 0
    return ((num_of_field == EOF ) ?
        ((int32_t) SLIP_INCORRECT_INPUT) : num_of_field) ;
}

/* Safely write the output as a null terminated string in a block of memory,
 * which is pointed to by a pointer stored in str. The block of memory must be
 * freed using mpfr_free_str. The return value is the number of characters
 * written in the string, excluding the null-terminator, or a negative value if
 * an error occurred */
 #if 0
/* This function is currently unused, but kept here for future reference. */
int32_t SLIP_mpfr_asprintf (char **str, const char *template, ... )
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpfr_vasprintf
    int32_t num_of_char;
    va_list args;
    va_start(args, template);
    num_of_char = mpfr_vasprintf(str, template, args);
    va_end (args);

    // Finish the wrapper and return num_of_char if successful
    SLIP_GMP_WRAPPER_FINISH;
    // mpfr_vasprintf returns a negative value if an error occurred
    if (num_of_char < 0)
    {
        return (int32_t) SLIP_INCORRECT_INPUT;
    }
    else
    {
        return num_of_char;
    }
}

/* Safely free a string allocated. */
SLIP_info SLIP_mpfr_free_str (char *str)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpfr_free_str
    mpfr_free_str(str);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}
#endif

/* Safely print to the stream fp. Return positive value (the number of
 * characters written) upon success, otherwise return negative value (error
 * code) */
int32_t SLIP_mpfr_fprintf(FILE *fp, const char *format, ... )
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpfr_vfprintf
    int32_t num_of_char;
    va_list args;
    va_start(args, format);
    num_of_char = mpfr_vfprintf(fp, format, args);
    va_end (args);
    // Free cache from mpfr_vfprintf. Even though mpfr_free_cache is
    // called in SLIP_LU_final(), it has to be called here to
    // prevent memory leak in some rare situations.
    mpfr_free_cache();


    // Finish the wrapper and return num_of_char if successful
    SLIP_GMP_WRAPPER_FINISH;
    // mpfr_vfprintf returns -1 if an error occurred.
    return ((num_of_char < 0) ?
        ((int32_t) SLIP_INCORRECT_INPUT) : num_of_char) ;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Integer (mpz_t type) functions-----------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


/* Purpose: Safely initialize an mpz_t number */
SLIP_info SLIP_mpz_init
(
    mpz_t x
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_init
    mpz_init(x);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}


/* Purpose: Safely initialize an mpz_t number with space for size bits */
SLIP_info SLIP_mpz_init2
(
    mpz_t x,                // Number to be initialized
    const uint64_t size     // size of the number
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_init2
    mpz_init2(x, (unsigned long int) size);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpz number = to an mpz number, i.e., x = y */
SLIP_info SLIP_mpz_set
(
    mpz_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_set
    mpz_set(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpz number = to an unsigned int, i.e., x = y */
SLIP_info SLIP_mpz_set_ui
(
    mpz_t x,
    const uint64_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_set_ui
    mpz_set_ui(x, (unsigned long int) y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpz number = a signed int */
SLIP_info SLIP_mpz_set_si
(
    mpz_t x,
    const int32_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_set_si
    mpz_set_si(x, (int) y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpz number = a double */
#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info SLIP_mpz_set_d
(
    mpz_t x,
    const double y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_set_d
    mpz_set_d(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set a double number = a mpz */
SLIP_info SLIP_mpz_get_d
(
    double *x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpz_set_d
    *x = mpz_get_d(y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}
#endif

/* Purpose: Safely set an mpz number = mpq number */
SLIP_info SLIP_mpz_set_q
(
    mpz_t x,
    const mpq_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_set_q
    mpz_set_q(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compute a = b*c */
SLIP_info SLIP_mpz_mul
(
    mpz_t a,
    const mpz_t b,
    const mpz_t c
)
{
    // Start the GMP Wrapper
    SLIP_GMPZ_WRAPPER_START(a);

    // Call gmp mul
    mpz_mul(a, b, c);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compute a = b+c */
#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info SLIP_mpz_add
(
    mpz_t a,
    const mpz_t b,
    const mpz_t c
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(a);

    // call mpz_add
    mpz_add(a,b,c);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpz number += product of two mpz numbers,
 * i.e., x = x + y*z */
/* This function is currently unused, but kept here for future reference. */
SLIP_info SLIP_mpz_addmul
(
    mpz_t x,
    const mpz_t y,
    const mpz_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_addmul
    mpz_addmul(x, y, z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}
#endif

/* Purpose: Safely set an mpz number = itself minus a product of
 * mpz numbers, i.e., x = x - y*z
 */
SLIP_info SLIP_mpz_submul
(
    mpz_t x,
    const mpz_t y,
    const mpz_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_submul
    mpz_submul(x, y, z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safe version of exact integer division, i.e., x = y / z */
SLIP_info SLIP_mpz_divexact
(
    mpz_t x,
    const mpz_t y,
    const mpz_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_divexact
    mpz_divexact(x, y, z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compute the gcd of two mpz_t numbers, i.e., x = gcd(y, z) */
SLIP_info SLIP_mpz_gcd
(
    mpz_t x,
    const mpz_t y,
    const mpz_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_gcd
    mpz_gcd(x, y, z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compute the lcm of two mpz numbers */
SLIP_info SLIP_mpz_lcm
(
    mpz_t lcm,   // lcm of x and y
    const mpz_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(lcm);

    // call mpz_lcm
    mpz_lcm(lcm, x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set x = |y| */
SLIP_info SLIP_mpz_abs
(
    mpz_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpz_abs
    mpz_abs(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compare two mpz numbers,
 * r > 0 if x > y, r = 0 if x = y, and r < 0 if x < y */
SLIP_info SLIP_mpz_cmp
(
    int32_t *r,
    const mpz_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpz_cmp
    *r = (int32_t) mpz_cmp(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compare the absolute value of two mpz numbers,
 * r > 0 if |x| > |y|, r = 0 if |x| = |y|, and r < 0 if |x| < |y| */
SLIP_info SLIP_mpz_cmpabs
(
    int32_t *r,
    const mpz_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpz_cmpabs
    *r = (int32_t) mpz_cmpabs(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compare a mpz number with a uint64_t integer
 * r > 0 if x > y, r = 0 if x = y, and r < 0 if x < y */
SLIP_info SLIP_mpz_cmp_ui
(
    int32_t *r,
    const mpz_t x,
    const uint64_t y
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpz_cmp_ui
    *r = (int32_t) mpz_cmp_ui(x, (unsigned long int) y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set sgn = 0 if x = 0, otherwise, sgn = x/|x| */
SLIP_info SLIP_mpz_sgn
(
    int32_t *sgn,
    const mpz_t x
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpz_sgn
    *sgn = (int32_t) mpz_sgn(x);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely return the size of x measured in number of digits
 * in the given base */
SLIP_info SLIP_mpz_sizeinbase
(
    uint64_t *size,
    const mpz_t x,
    int32_t base
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpz_lcm
    *size = (uint64_t) mpz_sizeinbase(x, (int) base);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Rational (mpq type) functions------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/* Purpose: Safely initialize an mpq_t number */
SLIP_info SLIP_mpq_init
(
    mpq_t x
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_init
    mpq_init(x);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpq number = to an mpq number, i.e., x = y */
SLIP_info SLIP_mpq_set
(
    mpq_t x,
    const mpq_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_set
    mpq_set(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpq number = an mpz number. i.e., x = y */
SLIP_info SLIP_mpq_set_z
(
    mpq_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_set_z
    mpq_set_z(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpq number = a double */
SLIP_info SLIP_mpq_set_d
(
    mpq_t x,
    const double y
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_set_d
    mpq_set_d(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpq number as the fraction of two
 * unsigned ints. i.e., x = y / z
 */
SLIP_info SLIP_mpq_set_ui
(
    mpq_t x,
    const uint64_t y,
    const uint64_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_set_ui
    mpq_set_ui(x, (unsigned long int) y, (unsigned long int) z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set the numerator of an mpq number */
SLIP_info SLIP_mpq_set_num
(
    mpq_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_set_num
    mpq_set_num(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set the denominator of an mpq number */
SLIP_info SLIP_mpq_set_den
(
    mpq_t x,
    const mpz_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_set_den
    mpq_set_den(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpz number = denominator of an mpq number */
SLIP_info SLIP_mpq_get_den
(
    mpz_t x,
    const mpq_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpq_get_den
    mpq_get_den(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set a double = a mpq number*/
SLIP_info SLIP_mpq_get_d
(
    double *x,
    const mpq_t y
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpq_get_d
    *x = mpq_get_d(y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpq number = absolute value of mpq */
SLIP_info SLIP_mpq_abs
(
    mpq_t x,
    const mpq_t y
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_abs
    mpq_abs(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely add two mpq numbers, i.e., x = y+z */
SLIP_info SLIP_mpq_add
(
    mpq_t x,
    const mpq_t y,
    const mpq_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_add
    mpq_add(x, y, z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely multiply two mpq numbers, i.e., x = y*z */
SLIP_info SLIP_mpq_mul
(
    mpq_t x,
    const mpq_t y,
    const mpq_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_mul
    mpq_mul(x, y, z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely divide two mpq numbers, i.e., x = y/z */
SLIP_info SLIP_mpq_div
(
    mpq_t x,
    const mpq_t y,
    const mpq_t z
)
{
    // Start the GMP wrapper
    SLIP_GMPQ_WRAPPER_START(x);

    // call mpq_div
    mpq_div(x, y, z);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compare two mpq numbers,
 * r > 0 if x > y, r = 0 if x = y, and r < 0 if x < y */
SLIP_info SLIP_mpq_cmp
(
    int32_t *r,
    const mpq_t x,
    const mpq_t y
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpq_cmp
    *r = (int32_t) mpq_cmp(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely compare x and num/den. r > 0 if x > num/den,
 * r = 0 if x = num/den, and r < 0 if x < num/den */
SLIP_info SLIP_mpq_cmp_ui
(
    int32_t *r,
    const mpq_t x,
    const uint64_t num,
    const uint64_t den
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpq_cmp_ui
    *r = (int32_t) mpq_cmp_ui(x, (unsigned long int) num,
                                 (unsigned long int) den);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely check if two mpq numbers equal,
 * r = 0 (r = false) if x != y, r != 0 (r = true) if x = y */
SLIP_info SLIP_mpq_equal
(
    int32_t *r,
    const mpq_t x,
    const mpq_t y
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpq_equal
    *r = (int32_t) mpq_equal(x, y);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Floating Point (mpfr type) functions-----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/* Purpose: Safely initialize an mpfr_t number */
SLIP_info SLIP_mpfr_init2
(
    mpfr_t x,       // Floating point number to initialize
    uint64_t size    // # of bits in x
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_init2
    mpfr_init2(x, (unsigned long int) size);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpfr number = to an mpfr number, i.e., x = y */
#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info SLIP_mpfr_set
(
    mpfr_t x,
    const mpfr_t y,
    const mpfr_rnd_t rnd
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_set
    mpfr_set(x, y, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}
#endif

/* Purpose: Safely set an mpfr number = to a double, i.e., x = y */
SLIP_info SLIP_mpfr_set_d
(
    mpfr_t x,
    const double y,
    const mpfr_rnd_t rnd  // MPFR rounding scheme used
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_set_d
    mpfr_set_d(x, y, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpfr number = to an mpq number */
SLIP_info SLIP_mpfr_set_q
(
    mpfr_t x,
    const mpq_t y,
    const mpfr_rnd_t rnd
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_set_q
    mpfr_set_q(x, y, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpfr number = to an mpz number */
SLIP_info SLIP_mpfr_set_z
(
    mpfr_t x,
    const mpz_t y,
    const mpfr_rnd_t rnd
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_set_q
    mpfr_set_z(x, y, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpz number = to an mpfr number, i.e., x = y */
SLIP_info SLIP_mpfr_get_z
(
    mpz_t x,
    const mpfr_t y,
    const mpfr_rnd_t rnd  // MPFR rounding scheme used
)
{
    // Start the GMP wrapper
    SLIP_GMPZ_WRAPPER_START(x);

    // call mpfr_get_z
    mpfr_get_z(x, y, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set a double = to a mpfr number, i.e., x = y */
SLIP_info SLIP_mpfr_get_d
(
    double *x,
    const mpfr_t y,
    const mpfr_rnd_t rnd  // MPFR rounding scheme used
)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpfr_get_d
    *x = mpfr_get_d(y, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely multiply mpfr numbers, x = y*z */
SLIP_info SLIP_mpfr_mul
(
    mpfr_t x,
    const mpfr_t y,
    const mpfr_t z,
    const mpfr_rnd_t rnd  // MPFR rounding mode
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_mul
    mpfr_mul(x, y, z, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpfr number = to a product of an mpfr_t and double,
 * i.e., x = y*z
 */
SLIP_info SLIP_mpfr_mul_d
(
    mpfr_t x,
    const mpfr_t y,
    const double z,
    const mpfr_rnd_t rnd  // MPFR rounding scheme used
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_mul_d
    mpfr_mul_d(x, y, z, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set a mpfr number = a mpfr number divided by a double,
 * i.e., x = y/z
 */
SLIP_info SLIP_mpfr_div_d
(
    mpfr_t x,
    const mpfr_t y,
    const double z,
    const mpfr_rnd_t rnd  // MPFR rounding scheme used
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_mul_d
    mpfr_div_d(x, y, z, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely set an mpfr number = power of two ints, i.e.,
 * x = y^z
 */
SLIP_info SLIP_mpfr_ui_pow_ui
(
    mpfr_t x,
    const uint64_t y,
    const uint64_t z,
    const mpfr_rnd_t rnd  // MPFR rounding mode
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpfr_ui_pow_ui
    mpfr_ui_pow_ui(x, (unsigned long int) y, (unsigned long int) z, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely take the log2 of an mpfr number */
SLIP_info SLIP_mpfr_log2
(
    mpfr_t x,
    const mpfr_t y,
    const mpfr_rnd_t rnd
)
{
    // Start the GMP wrapper
    SLIP_GMPFR_WRAPPER_START(x);

    // call mpq_add
    mpfr_log2(x, y, rnd);

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

/* Purpose: Safely free all caches and pools used by MPFR internally */
SLIP_info SLIP_mpfr_free_cache(void)
{
    // Start the GMP wrapper
    SLIP_GMP_WRAPPER_START;

    // call mpfr_free_cache
    mpfr_free_cache();

    // Finish the wrapper and return 0 if successful
    SLIP_GMP_WRAPPER_FINISH;
    return SLIP_OK;
}

