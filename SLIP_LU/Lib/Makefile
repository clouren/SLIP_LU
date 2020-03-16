#-------------------------------------------------------------------------------
# SLIP_LU/Lib/Makefile
#-------------------------------------------------------------------------------

# To just compile the library, do:
#       make
#
# To run a demo using the library
#       cd ../Demo ; make

LIBRARY = libsliplu
VERSION = 1.0.0
SO_VERSION = 1

default: library

include ../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
# CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SLIP_LU depends on SuiteSparse_config, AMD, COLAMD, M, GMP, and MPFR
LDLIBS += -lsuitesparseconfig -lamd -lcolamd -lm -lgmp -lmpfr

C = $(CC) $(CF) -I../Include -I../../COLAMD/Include -I../../AMD/Include -I../../SuiteSparse_config

all: install

OBJ = \
    slip_array_div.o slip_array_mul.o slip_back_sub.o \
    SLIP_build_dense_mpz.o \
    SLIP_build_dense_mpq.o \
    SLIP_build_dense_mpfr.o  \
    SLIP_build_dense_int.o \
    SLIP_build_dense_double.o \
    SLIP_build_sparse_trip_mpz.o \
    SLIP_build_sparse_trip_mpq.o \
    SLIP_build_sparse_trip_mpfr.o \
    SLIP_build_sparse_trip_int.o \
    SLIP_build_sparse_trip_double.o \
    SLIP_build_sparse_csc_mpz.o \
    SLIP_build_sparse_csc_mpq.o \
    SLIP_build_sparse_csc_mpfr.o \
    SLIP_build_sparse_csc_double.o \
    SLIP_build_sparse_csc_int.o \
    SLIP_calloc.o SLIP_check_solution.o slip_create_LU_analysis.o \
    slip_create_dense.o SLIP_create_double_mat.o SLIP_delete_int_mat.o \
    SLIP_create_int_mat.o SLIP_create_mpfr_array.o SLIP_create_mpfr_mat.o \
    SLIP_create_mpq_array.o SLIP_create_mpq_mat.o slip_create_mpz_array2.o \
    SLIP_create_mpz_array.o SLIP_create_mpz_mat.o \
    slip_create_sparse.o slip_cumsum.o SLIP_delete_dense.o \
    SLIP_delete_double_mat.o SLIP_delete_LU_analysis.o \
    SLIP_delete_mpfr_array.o SLIP_solve_mpfr.o SLIP_create_default_options.o\
    SLIP_delete_mpfr_mat.o SLIP_delete_mpq_array.o SLIP_delete_mpq_mat.o \
    SLIP_delete_mpz_array.o SLIP_delete_mpz_mat.o SLIP_delete_sparse.o \
    slip_dense_alloc.o slip_dfs.o slip_expand_double_array.o \
    slip_expand_double_mat.o slip_expand_mpfr_array.o slip_expand_mpfr_mat.o \
    slip_expand_mpq_array.o slip_expand_mpq_mat.o SLIP_finalize.o \
    slip_forward_sub.o SLIP_free.o slip_get_column.o SLIP_get_double_soln.o \
    slip_get_largest_pivot.o SLIP_get_mpfr_soln.o slip_get_nonzero_pivot.o \
    slip_get_pivot.o slip_get_smallest_pivot.o SLIP_gmp.o SLIP_initialize.o \
    SLIP_initialize_expert.o SLIP_LU_analyze.o SLIP_LU_factorize.o \
    slip_lu_info.o \
    SLIP_LU_solve.o SLIP_malloc.o slip_mpz_populate_mat.o SLIP_permute_x.o \
    slip_reach.o SLIP_realloc.o slip_ref_triangular_solve.o \
    slip_reset_int_array2.o \
    slip_reset_int_array.o slip_reset_mpz_array.o SLIP_scale_x.o \
    SLIP_solve_double.o \
    SLIP_solve_mpq.o slip_sort_xi.o slip_sparse_alloc2.o slip_sparse_alloc.o \
    slip_sparse_collapse.o slip_sparse_realloc.o SLIP_spok.o \
    slip_trip_to_mat.o

$(OBJ): ../Include/SLIP_LU.h ../Source/SLIP_LU_internal.h Makefile

%.o: ../Source/%.c ../Include/SLIP_LU.h ../Source/SLIP_LU_internal.h
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------
# Create the din
#-------------------------------------------------------------------------------

library:
	$(MAKE) install INSTALL=$(SUITESPARSE)

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $(AR_TARGET)

#-------------------------------------------------------------------------------
# install (shared C version only)
#-------------------------------------------------------------------------------

# install SLIP_LUI
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CC) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/SLIP_LU.h $(INSTALL_INCLUDE)
	$(CP) ../Doc/SLIP_LU_UserGuide.pdf $(INSTALL_DOC)
	$(CP) ../README.md $(INSTALL_DOC)/SLIP_LU_README.md
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/SLIP_LU.h
	chmod 644 $(INSTALL_DOC)/SLIP_LU_UserGuide.pdf
	chmod 644 $(INSTALL_DOC)/SLIP_LU_README.md

# uninstall SLIP_LU
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SLIP_LU.h
	$(RM) $(INSTALL_DOC)/SLIP_LU_UserGuide.pdf
	$(RM) $(INSTALL_DOC)/SLIP_LU_README.md

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) -r $(CLEAN)

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

