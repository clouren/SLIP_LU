#-------------------------------------------------------------------
# Compiler selection and flags
#--------------------------------------------------------------------
CC  = gcc
CCC = g++
debug = -g -o
optimize = -O2 -o


#--------------------------------------------------------------------
#Library flags
#--------------------------------------------------------------------
CCLSTD = -std=c++11
CCLGMP =   -lgmp -lgmpxx
CCLMPFR = -lmpfr
CCLNCOLS = -lcolamd -lamd
#
#

#--------------------------------------------------------------------
# Applications
#--------------------------------------------------------------------- 
APPS = SLIP_LU 
Debugs = debugSLIP

all: $(APPS) $(Debugs)
debug: $(Debugs)
	
SLIP_LU: SLIP_LU.cpp
	$(CCC) SLIP_LU.cpp $(optimize) SLIP_LU.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) 

debugSLIP: SLIP_LU.cpp
	$(CCC) SLIP_LU.cpp $(debug) debug_SLIP_LU.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD)

	

#--------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------
clean:
	rm *.exe
