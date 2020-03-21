#-------------------------------------------------------------------------------
# SLIP_LU/Tcov/Makefile: compile and run SLIP_LU test coverage
#-------------------------------------------------------------------------------

# SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
# Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
# SLIP_LU/License for the license.

#-------------------------------------------------------------------------------

default: run

include ../../SuiteSparse_config/SuiteSparse_config.mk

# Linux test coverage
CC = gcc
CFLAGS = -g -fprofile-arcs -ftest-coverage \
        -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
        -Wredundant-decls -Wnested-externs -Wdisabled-optimization -std=c99 \
        -Wno-unused-parameter -I../Include -I../Source -I../Demo \
        -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include \
        -DSLIP_GMP_LIST_INIT=2 -DSLIP_LU_TCOV

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd

# run all statement coverage tests, and then check for 100% coverage
run: runtests
	./covall

all: test SLIPLU

CS = \
    slip_array_div.o slip_array_mul.o slip_back_sub.o \
    SLIP_build_dense_mpz.o \
    SLIP_build_dense_mpq.o \
    SLIP_build_dense_mpfr.o \
    SLIP_build_dense_int32.o \
    SLIP_build_dense_double.o \
    SLIP_build_sparse_trip_mpz.o \
    SLIP_build_sparse_trip_mpq.o \
    SLIP_build_sparse_trip_mpfr.o \
    SLIP_build_sparse_trip_int32.o \
    SLIP_build_sparse_trip_double.o \
    SLIP_build_sparse_csc_mpz.o \
    SLIP_build_sparse_csc_mpq.o \
    SLIP_build_sparse_csc_mpfr.o  \
    SLIP_build_sparse_csc_int32.o \
    SLIP_build_sparse_csc_double.o \
    SLIP_check_solution.o slip_create_LU_analysis.o slip_create_dense.o \
    SLIP_create_double_mat.o SLIP_delete_int32_mat.o SLIP_free.o \
    SLIP_create_int32_mat.o SLIP_create_mpfr_array.o SLIP_create_mpfr_mat.o \
    SLIP_create_mpq_array.o SLIP_create_mpq_mat.o slip_create_mpz_array2.o \
    SLIP_create_mpz_array.o SLIP_create_mpz_mat.o \
    slip_create_sparse.o slip_cumsum.o SLIP_delete_dense.o \
    SLIP_delete_double_mat.o SLIP_delete_LU_analysis.o \
    SLIP_delete_mpfr_array.o SLIP_solve_mpfr.o  SLIP_create_default_options.o\
    SLIP_delete_mpfr_mat.o SLIP_delete_mpq_array.o SLIP_delete_mpq_mat.o \
    SLIP_delete_mpz_array.o SLIP_delete_mpz_mat.o SLIP_delete_sparse.o \
    slip_dense_alloc.o slip_dfs.o slip_expand_double_array.o \
    slip_expand_double_mat.o slip_expand_mpfr_array.o slip_expand_mpfr_mat.o \
    slip_expand_mpq_array.o slip_expand_mpq_mat.o SLIP_finalize.o \
    slip_forward_sub.o slip_get_column.o SLIP_get_double_soln.o \
    slip_get_largest_pivot.o SLIP_get_mpfr_soln.o slip_get_nonzero_pivot.o \
    slip_get_pivot.o slip_get_smallest_pivot.o SLIP_gmp.o SLIP_initialize.o \
    SLIP_LU_analyze.o SLIP_LU_factorize.o slip_lu_info.o \
    SLIP_LU_solve.o slip_ref_triangular_solve.o \
    slip_mpz_populate_mat.o SLIP_permute_x.o slip_reach.o \
    slip_reset_int32_array2.o slip_reset_int32_array.o \
    slip_reset_mpz_array.o SLIP_scale_x.o SLIP_solve_double.o \
    SLIP_solve_mpq.o slip_sort_xi.o slip_sparse_alloc2.o slip_sparse_alloc.o \
    slip_sparse_collapse.o slip_sparse_realloc.o SLIP_spok.o \
    slip_trip_to_mat.o tcov_malloc_test.o SLIP_initialize_expert.o

$(CS): ../Include/SLIP_LU.h ../Source/SLIP_LU_internal.h tcov_malloc_test.h ../Demo/demos.h

.PRECIOUS: SLIP_%.c slip_%.c SLIPLU.c demos.c

SLIP_%.c:
	- ln -s ../Source/$@

slip_%.c:
	- ln -s ../Source/$@

SLIPLU.c:
	- ln -s ../Demo/$@

demos.c:
	- ln -s ../Demo/$@

test: $(CS) cov_test.c demos.c
	$(CC) $(LDFLAGS) cov_test.c demos.c $(CFLAGS) -o cov_test $(CS) $(LDLIBS)

SLIPLU: $(CS) SLIPLU.c demos.c
	$(CC) $(LDFLAGS) SLIPLU.c demos.c $(CFLAGS) -o SLIPLU $(CS) $(LDLIBS)

# run all statement coverage tests
runtests: all
	- ./SLIPLU p 2 q 0
	- ./SLIPLU p 3 q 1 o2 1
	- ./SLIPLU p 4 q 2 o2 1 f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- ./SLIPLU p 5
	- ./cov_test
	- ./cov_test 0 1 1 1

# To run with valgrind:
V = valgrind

# run all statement coverage tests but with valgrind
vtests: all
	- $(V) ./SLIPLU p 2 q 0
	- $(V) ./SLIPLU p 3 q 1 o2 1
	- $(V) ./SLIPLU p 4 q 2 o2 1 f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- $(V) ./SLIPLU p 5
	- $(V) ./cov_test
	- $(V) ./cov_test 0 1 1 1

# remove all files not in the original distribution, including symbolic links
clean:
	- $(RM) *.o *.bbg *.da *.gcov *.gcda *gcno
	- $(RM) cov_test SLIPLU *.out *.a cov.sort out
	- $(RM) -r SLIP_*.c slip_*.c *.dSYM SLIPLU.c demos.c

purge: distclean

distclean: clean

