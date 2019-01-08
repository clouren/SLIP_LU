This software package is used to solve a sparse systems of linear equations exactly using the SLIP LU factorization.


This folder contains eight C++ files: SLIP_LU, SLIP_random, example, example2, example3, example4, example5
This folder also contains header files (including Suitesparse headers). Do not modify these headers.
This folder also contains MATLAB interface

*********SLIP_LU*********
Purpose: Exactly solve a sparse system of linear equations using a given input matrix and right hand side vector file. This code can 
         output the final solution to a user specified output file in either double precision or full precision rational numbers.
         If you intend to use SLIP LU within another program, please refer to examples for help with this.

USAGE: ./SLIP_LU.exe Followed by:
    f or file: Filename.         SLIP LU expects the format to be f MATRIX_NAME RHS_NAME. The matrix must be stored in matrix market format 
                                 Please refer to http://math.nist.gov/MatrixMarket/formats.html for information on matrix market format. The 
                                 right hand side vector must be stored as a dense vector.
            
    of or outfile: Outfile name: SLIP LU expects the formalt to be of OUTNAME. 
                             
    c or check: Check param:     This indicates that SLIP LU will check the solution of the problem. 
                                 WARNING: This is slow and should only be used for verification

    p or piv: Pivot param:       This parameter tells SLIP LU what type of pivoting to use. The options are as follows:
                                 0: Smallest pivot: Default and recommended
                                 1: Diagonal pivoting
                                 2: First nonzero per column chosen as pivot
                                 3: Diagonal pivoting with tolerance for smallest pivot
                                 4: Diagonal pivoting with tolerance for largest pivot
                                 5: Largest pivot
                                 
    q or col: Column order param:This tells SLIP LU what column ordering to use. Options are:
                                 0: COLAMD: Default
                                 1: AMD
                                 2: None: Not recommended for sparse matrices
                                 3: UMFPACK P and Q
    
    
    t or tol: tolerance param:   This is SLIP LU's tolerance. Only necessary if some sort of tolerance pivoting is used
    
    out2 or o2: output param:    This indicates SLIP LU will output information about the ordering, SLIP LU
    
    out or o: output param:      This indicates SLIP LU will output the solution to a file. This should be followed by one of the following integers
                                 1: The solution will be output in full precision rational arithmetic
                                 2: The solution will be output in double precision 
                                 3: The solution will be output in some user specified precision. This must be input as o2 3 USER_PRECISION. Precision must be an int
                                 
                                 
***********NOTE: If none of these are used, SLIP LU uses default parameters defined in SLIP_LU_config.h********

*********SLIP_random*********
Purpose: Compare SLIP_LU with REF_LU on dense, randomly generated matrices. There are a different set of command line arguments for this one.
USAGE: ./SLIP_random.exe Followed by:
	n:             dimension of matrix to generate
	
	b  or numRHS:  number of right hand side vectors to generate
	
	l  or lower:   lower bound of randomly generated numbers
	
	u  or upper:   upper bound of entries to generate
	
	s1 or seed1:   first random number seed used to generate entries
	
	s2 or seed2:   second random number seed used to generate entries
	
	d  or density: density of matrix generated


*********example*********
Purpose: Demonstrate the use of SLIP LU for a user given input matrix

*********example2*********
Purpose: Demonstrate the use of SLIP LU for a matrix to be read in

*********example3*********
Purpose: Demonstrate the use of SLIP LU for a matrix stored in mpfr precision

*********example4*********
Purpose: Demonstrate the use of SLIP LU for multiple RHS vectors

*********example5*********
Purpose: Demonstrate the use of SLIP LU for multiple RHS vectors read in from file
