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
CCLNUMF = -lumfpack
#
#

#--------------------------------------------------------------------
# Applications
#--------------------------------------------------------------------- 
APPS = SLIP_LU SLIP_LU_random SLIP_chol SLIP_mg example example2 example3 example4 example5
LU = SLIP_LU
Debugs = debugSLIP debugChol debugRandom debugexample debugexample2 debugexample3 debugexample4 debugexample5 debugMG

all: $(APPS) $(Debugs)
lu: $(LUS)
debug: $(Debugs)
run: $(tests)
	
SLIP_LU: SLIP_LU.cpp
	$(CCC) SLIP_LU.cpp $(optimize) SLIP_LU.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

debugSLIP: SLIP_LU.cpp
	$(CCC) SLIP_LU.cpp $(debug) debug_SLIP_LU.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

SLIP_chol: SLIP_chol.cpp
	$(CCC) SLIP_chol.cpp $(optimize) SLIP_chol.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
debugChol: SLIP_chol.cpp
	$(CCC) SLIP_chol.cpp $(debug) debug_SLIP_chol.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
SLIP_mg: SLIP_mg.cpp
	$(CCC) SLIP_mg.cpp $(optimize) SLIP_mg.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

debugMG: SLIP_mg.cpp
	$(CCC) SLIP_mg.cpp $(debug) debug_SLIP_mg.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
SLIP_LU_random: SLIP_LU_random.cpp
	$(CCC) SLIP_LU_random.cpp $(optimize) SLIP_random.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
debugRandom: SLIP_LU_random.cpp
	$(CCC) SLIP_LU_random.cpp $(debug) debug_random.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)

example: ./Demos/example.cpp
	$(CCC) ./Demos/example.cpp $(optimize) ./Demos/example.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
debugexample: ./Demos/example.cpp
	$(CCC) ./Demos/example.cpp $(debug) ./Demos/debug_example.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
example2: ./Demos/example2.cpp
	$(CCC) ./Demos/example2.cpp $(optimize) ./Demos/example2.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
debugexample2: ./Demos/example2.cpp
	$(CCC) ./Demos/example2.cpp $(debug) ./Demos/debug_example2.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
example3: ./Demos/example3.cpp
	$(CCC) ./Demos/example3.cpp $(optimize) ./Demos/example3.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF) 
	
debugexample3: ./Demos/example3.cpp
	$(CCC) ./Demos/example3.cpp $(debug) ./Demos/debug_example3.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
example4: ./Demos/example4.cpp
	$(CCC) ./Demos/example4.cpp $(optimize) ./Demos/example4.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
debugexample4: ./Demos/example4.cpp
	$(CCC) ./Demos/example4.cpp $(debug) ./Demos/debug_example4.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
example5: ./Demos/example5.cpp
	$(CCC) ./Demos/example5.cpp $(optimize) ./Demos/example5.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	
debugexample5: ./Demos/example5.cpp
	$(CCC) ./Demos/example5.cpp $(debug) ./Demos/debug_example5.exe $(CCLGMP) $(CCLMPFR) $(CCLNCOLS) $(CCLSTD) $(CCLNUMF)
	

#--------------------------------------------------------------------
# Efficient testing after make
#---------------------------------------------------------------------
	
test: 
	./SLIP_LU.exe c; ./SLIP_LU.exe f ../../BasisLIB_INT/NSR8K.mat ../../BasisLIB_INT/NSR8K.v c; 
	cd Demos; ./example.exe; ./example2.exe; ./example3.exe; ./example4.exe; ./example5.exe; cd ..;
	./SLIP_random.exe; echo "All examples successful";


#--------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------
clean:
	rm *.exe
