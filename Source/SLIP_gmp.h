#ifndef SLIP_gmp
#define SLIP_gmp

#include <stdarg.h>
#include "SLIP_LU_internal.h"

// The GMP library does not handle out-of-memory failures.  However, it does
// provide a mechanism for passing function pointers that replace GMP's use of
// malloc, realloc, and free.  This mechanism is used to provide a try/catch
// mechanism for memory allocation errors, using setjmp and longjmp.

// When a GMP function is called, this wrapper keeps track of a list of objects
// allocated by that function.  The list is started fresh each time a GMP
// function is called.  If any allocation fails, the NULL pointer is not
// returned to GMP.  Instead, all allocated blocks in the list are freed,
// and slip_gmp_allocate returns directly to wrapper.


#ifndef SLIP_GMP_LIST_INIT
// A size of 32 ensures that the list never needs to be increased in size.
// The test coverage suite in SLIP_LU/Tcov reduces this initial size to
// exercise the code, in SLIP_LU/Tcov/Makefile.
#define SLIP_GMP_LIST_INIT 32
#endif

// uncomment this to print memory debugging info
// #define SLIP_GMP_MEMORY_DEBUG

#ifdef SLIP_GMP_MEMORY_DEBUG
void slip_gmp_dump ( void ) ;
#endif


#define SLIP_GMP_WRAPPER_START                                          \
{                                                                       \
    slip_gmp_nmalloc = 0 ;                                              \
    /* setjmp returns 0 if called from here, or > 0 if from longjmp */  \
    int32_t slip_gmp_status = setjmp (slip_gmp_environment) ;           \
    if (slip_gmp_status != 0)                                           \
    {                                                                   \
        /* failure from longjmp */                                      \
        slip_gmp_failure (slip_gmp_status) ;                            \
        return (SLIP_OUT_OF_MEMORY) ;                                   \
    }                                                                   \
}
#define SLIP_GMPZ_WRAPPER_START(x)                                      \
{                                                                       \
    slip_gmpz_archive = (mpz_t *) x;                                    \
    slip_gmpq_archive = NULL;                                           \
    slip_gmpfr_archive = NULL;                                          \
    SLIP_GMP_WRAPPER_START;                                             \
}

#define SLIP_GMPQ_WRAPPER_START(x)                                      \
{                                                                       \
    slip_gmpz_archive = NULL;                                           \
    slip_gmpq_archive =(mpq_t *) x;                                     \
    slip_gmpfr_archive = NULL;                                          \
    SLIP_GMP_WRAPPER_START;                                             \
}

#define SLIP_GMPFR_WRAPPER_START(x)                                     \
{                                                                       \
    slip_gmpz_archive = NULL;                                           \
    slip_gmpq_archive = NULL;                                           \
    slip_gmpfr_archive = (mpfr_t *) x;                                  \
    SLIP_GMP_WRAPPER_START;                                             \
}

#define SLIP_GMP_WRAPPER_FINISH                                         \
{                                                                       \
    /* clear (but do not free) the list.  The caller must ensure */     \
    /* the result is eventually freed. */                               \
    slip_gmpz_archive = NULL ;                                          \
    slip_gmpq_archive = NULL ;                                          \
    slip_gmpfr_archive = NULL ;                                         \
    slip_gmp_nmalloc = 0 ;                                              \
}

// free a block of memory, and also remove it from the archive if it's there
#define SLIP_SAFE_FREE(p)                                               \
{                                                                       \
    if (slip_gmpz_archive != NULL)                                      \
    {                                                                   \
        if (p == MPZ_PTR(*slip_gmpz_archive))                           \
        {                                                               \
            MPZ_PTR(*slip_gmpz_archive) = NULL ;                        \
        }                                                               \
    }                                                                   \
    else if (slip_gmpq_archive != NULL)                                 \
    {                                                                   \
        if (p == MPZ_PTR(MPQ_NUM(*slip_gmpq_archive)))                  \
        {                                                               \
            MPZ_PTR(MPQ_NUM(*slip_gmpq_archive)) = NULL ;               \
        }                                                               \
        if (p == MPZ_PTR(MPQ_DEN(*slip_gmpq_archive)))                  \
        {                                                               \
            MPZ_PTR(MPQ_DEN(*slip_gmpq_archive)) = NULL ;               \
        }                                                               \
    }                                                                   \
    else if (slip_gmpfr_archive != NULL)                                \
    {                                                                   \
        if (p == MPFR_REAL_PTR(*slip_gmpfr_archive))                    \
        {                                                               \
            MPFR_MANT(*slip_gmpfr_archive) = NULL ;                     \
        }                                                               \
    }                                                                   \
    SLIP_FREE (p) ;                                                     \
}

extern int64_t slip_gmp_ntrials ;

SLIP_info slip_gmp_fprintf(FILE *fp, const char *format, ... );

SLIP_info slip_gmp_printf(const char *format, ... );

SLIP_info slip_gmp_fscanf(FILE *fp, const char *format, ... );

#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info slip_mpfr_asprintf (char **str, const char *template, ... );

SLIP_info slip_mpfr_free_str (char *str);
#endif

SLIP_info slip_mpfr_fprintf(FILE *fp, const char *format, ... );

bool slip_gmp_init (void) ;

void slip_gmp_finalize (void) ;
 
void *slip_gmp_allocate (size_t size) ;

void slip_gmp_free (void *p, size_t size) ;

void *slip_gmp_reallocate (void *p_old, size_t old_size, size_t new_size );

void slip_gmp_failure (int32_t status) ;

SLIP_info slip_mpz_init(mpz_t x) ;

SLIP_info slip_mpz_init2(mpz_t x, const uint64_t size) ;

SLIP_info slip_mpz_set(mpz_t x, const mpz_t y);

SLIP_info slip_mpz_set_ui(mpz_t x, const uint64_t y) ;

SLIP_info slip_mpz_set_si(mpz_t x, const int32_t y);

#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info slip_mpz_set_d(mpz_t x, const double y);
#endif

SLIP_info slip_mpz_get_d(double *x, const mpz_t y) ;

SLIP_info slip_mpz_set_q(mpz_t x, const mpq_t y) ;

SLIP_info slip_mpz_mul(mpz_t a, const mpz_t b, const mpz_t c) ;

#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info slip_mpz_add(mpz_t a, const mpz_t b, const mpz_t c);

SLIP_info slip_mpz_addmul(mpz_t x, const mpz_t y, const mpz_t z) ;
#endif

SLIP_info slip_mpz_submul(mpz_t x, const mpz_t y, const mpz_t z);

SLIP_info slip_mpz_divexact(mpz_t x, const mpz_t y, const mpz_t z) ;

SLIP_info slip_mpz_gcd(mpz_t x, const mpz_t y, const mpz_t z) ;

SLIP_info slip_mpz_lcm(mpz_t lcm, const mpz_t x, const mpz_t y) ;

SLIP_info slip_mpz_abs(mpz_t x, const mpz_t y) ;

SLIP_info slip_mpz_cmp(int32_t *r, const mpz_t x, const mpz_t y) ;

SLIP_info slip_mpz_cmpabs(int32_t *r, const mpz_t x, const mpz_t y) ;

SLIP_info slip_mpz_cmp_ui(int32_t *r, const mpz_t x, const uint64_t y) ;

SLIP_info slip_mpz_sgn(int32_t *sgn, const mpz_t x) ;

SLIP_info slip_mpz_sizeinbase(uint64_t *size, const mpz_t x, int32_t base) ;

SLIP_info slip_mpq_init(mpq_t x) ;

SLIP_info slip_mpq_set(mpq_t x, const mpq_t y);

SLIP_info slip_mpq_set_z(mpq_t x, const mpz_t y);

SLIP_info slip_mpq_set_d(mpq_t x, const double y) ;

SLIP_info slip_mpq_set_ui(mpq_t x, const uint64_t y, const uint64_t z);

SLIP_info slip_mpq_set_num(mpq_t x, const mpz_t y) ;

SLIP_info slip_mpq_set_den(mpq_t x, const mpz_t y);

SLIP_info slip_mpq_get_den(mpz_t x, const mpq_t y) ;

SLIP_info slip_mpq_get_d(double *x, const mpq_t y) ;

SLIP_info slip_mpq_abs(mpq_t x, const mpq_t y) ;

SLIP_info slip_mpq_add(mpq_t x, const mpq_t y, const mpq_t z) ;

SLIP_info slip_mpq_mul(mpq_t x, const mpq_t y, const mpq_t z) ;

SLIP_info slip_mpq_div(mpq_t x, const mpq_t y, const mpq_t z) ;

SLIP_info slip_mpq_cmp(int32_t *r, const mpq_t x, const mpq_t y) ;

SLIP_info slip_mpq_cmp_ui(int32_t *r, const mpq_t x,
                    const uint64_t num, const uint64_t den) ;

SLIP_info slip_mpq_equal(int32_t *r, const mpq_t x, const mpq_t y) ;

SLIP_info slip_mpfr_init2(mpfr_t x, const uint64_t size ) ;

#if 0
/* This function is currently unused, but kept here for future reference. */
SLIP_info slip_mpfr_set(mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;
#endif

SLIP_info slip_mpfr_set_d(mpfr_t x, const double y, const mpfr_rnd_t rnd) ;

SLIP_info slip_mpfr_set_q(mpfr_t x, const mpq_t y, const mpfr_rnd_t rnd ) ;

SLIP_info slip_mpfr_set_z(mpfr_t x, const mpz_t y, const mpfr_rnd_t rnd ) ;

SLIP_info slip_mpfr_get_z(mpz_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SLIP_info slip_mpfr_get_d(double *x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SLIP_info slip_mpfr_mul(mpfr_t x, const mpfr_t y, const mpfr_t z,
                    const mpfr_rnd_t rnd) ;

SLIP_info slip_mpfr_mul_d(mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SLIP_info slip_mpfr_div_d(mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SLIP_info slip_mpfr_ui_pow_ui(mpfr_t x, const uint64_t y, const uint64_t z,
                    const mpfr_rnd_t rnd) ;

SLIP_info slip_mpfr_log2(mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd );

SLIP_info slip_mpfr_free_cache(void);


#endif
