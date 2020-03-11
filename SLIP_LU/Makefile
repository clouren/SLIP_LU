#------------------------------------------------------------------------------
# SLIP_LU Makefile
#------------------------------------------------------------------------------

# SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
# Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
# SLIP_LU/License for the license.

#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: C

include ../SuiteSparse_config/SuiteSparse_config.mk

demos: all

C:
	( cd Lib ; $(MAKE) )
	( cd Demo ; $(MAKE) )

all: C cov

library:
	( cd Lib ; $(MAKE) )

# compile the static libraries only
static:
	( cd Lib    ; $(MAKE) static )

cov:
	( cd Tcov ; $(MAKE) )

clean:
	( cd Lib ; $(MAKE) clean )
	( cd Demo ; $(MAKE) clean )
	( cd Tcov ; $(MAKE) clean )
	- ( cd MATLAB/Source ; $(RM) *.o )
	- ( cd MATLAB    ; $(RM) *.o )
	( cd Doc ; $(MAKE) clean )

purge:
	( cd Lib ; $(MAKE) purge )
	( cd Demo ; $(MAKE) purge )
	( cd Tcov ; $(MAKE) purge )
	- ( cd MATLAB/Source ; $(RM) *.o *.mex* )
	- ( cd MATLAB    ; $(RM) *.o *.mex* )
	( cd Doc ; $(MAKE) clean )

distclean: purge

install: library
	( cd Lib ; $(MAKE) install )

uninstall:
	( cd Lib ; $(MAKE) uninstall )

