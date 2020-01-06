# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
#
# To run a demo using the library 
#       cd ../Demo ; make
#


LIBRARY = libslip
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
        
I = -I../Include -I../Source -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)
CP = cp -f

AR_TARGET = $(LIBRARY).a

all: install

CS =slip_array_div.o slip_array_mul.o slip_back_sub.o SLIP_build_dense_mpz.o \
    SLIP_build_dense_mpq.o SLIP_build_dense_mpfr.o SLIP_build_dense_int.o \
    SLIP_build_dense_double.o SLIP_build_sparse_trip_mpz.o SLIP_build_sparse_trip_mpq.o \
    SLIP_build_sparse_trip_mpfr.o SLIP_build_sparse_trip_int.o SLIP_build_sparse_trip_double.o \
    SLIP_build_sparse_ccf_mpz.o SLIP_build_sparse_ccf_mpq.o SLIP_build_sparse_ccf_mpfr.o \
    SLIP_build_sparse_ccf_double.o SLIP_build_sparse_ccf_int.o SLIP_calloc.o \
    SLIP_check_solution.o SLIP_create_LU_analysis.o  \
    SLIP_create_dense.o SLIP_create_double_mat.o SLIP_delete_int_mat.o \
    SLIP_create_int_mat.o SLIP_create_mpfr_array.o SLIP_create_mpfr_mat.o \
    SLIP_create_mpq_array.o SLIP_create_mpq_mat.o slip_create_mpz_array2.o \
    SLIP_create_mpz_array.o SLIP_create_mpz_mat.o \
    SLIP_create_sparse.o slip_cumsum.o SLIP_delete_dense.o \
    SLIP_delete_double_mat.o SLIP_delete_LU_analysis.o \
    SLIP_delete_mpfr_array.o SLIP_solve_mpfr.o  SLIP_create_default_options.o\
    SLIP_delete_mpfr_mat.o SLIP_delete_mpq_array.o SLIP_delete_mpq_mat.o \
    SLIP_delete_mpz_array.o SLIP_delete_mpz_mat.o SLIP_delete_sparse.o \
    slip_dense_alloc.o slip_dfs.o slip_expand_double_array.o \
    slip_expand_double_mat.o slip_expand_mpfr_array.o slip_expand_mpfr_mat.o \
    slip_expand_mpq_array.o slip_expand_mpq_mat.o SLIP_finalize.o \
    slip_forward_sub.o SLIP_free.o slip_get_column.o SLIP_get_double_soln.o \
    slip_get_largest_pivot.o SLIP_get_mpfr_soln.o slip_get_nonzero_pivot.o \
    slip_get_pivot.o slip_get_smallest_pivot.o SLIP_gmp.o SLIP_initialize.o \
    SLIP_initialize_expert.o SLIP_LU_analyze.o SLIP_LU_factorize.o slip_lu_info.o \
    SLIP_LU_solve.o SLIP_malloc.o slip_mpz_populate_mat.o SLIP_permute_x.o \
    slip_reach.o SLIP_realloc.o   slip_REF_triangular_solve.o slip_reset_int_array2.o \
    slip_reset_int_array.o slip_reset_mpz_array.o SLIP_scale_x.o SLIP_solve_double.o \
    SLIP_solve_mpq.o slip_sort_xi.o slip_sparse_alloc2.o slip_sparse_alloc.o \
    slip_sparse_collapse.o slip_sparse_realloc.o SLIP_spok.o \
    slip_trip_to_mat.o 

$(CS): ../Include/SLIP_LU.h ../Source/SLIP_LU_internal.h  Makefile

%.o: ../Source/%.c ../Include/SLIP_LU.h ../Source/SLIP_LU_internal.h
	$(CC) $(CF) $(I) -c $<

static: $(AR_TARGET)

$(AR_TARGET): $(CS)
	$(ARCHIVE)  $@ $^
	- $(RANLIB) $@

# install SLIP LU in this directory
install: $(AR_TARGET)

# uninstall SLIP LU
uninstall: purge

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) *.a *.obj *.dll *.dylib *.so *.so.*

