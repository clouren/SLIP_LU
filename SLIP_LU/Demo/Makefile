default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

# CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../Include -I../Source -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include

# LDFLAGS = -L../../lib

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd
CS = ../Lib/libslip.a $(LDLIBS)


all: lib example example2 example3 example4 example5 lu
	- ./example
	- ./example2
	- ./example3
	- ./example4
	- ./example5
	- ./SLIPLU

lib:
	( cd ../Lib ; $(MAKE) )

example: lib example.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example example.c $(CS)

example2: lib example2.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example2 example2.c demos.c $(CS)

example3: lib example3.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example3 example3.c  $(CS)

example4: lib example4.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example4 example4.c  $(CS)

example5: lib example5.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example5 example5.c demos.c $(CS)

lu: lib SLIPLU.c demos.h demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o SLIPLU SLIPLU.c demos.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example example2 example3 example4 example5 SLIPLU *.a *.dSYM *.obj *.dll
