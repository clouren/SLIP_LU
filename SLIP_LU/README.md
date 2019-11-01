This software package is used to solve a sparse systems of linear equations
exactly using the SLIP LU factorization. 


*********SLIPLU*********
Purpose: Exactly solve a sparse system of linear equations using a given input
         matrix and right hand side vector file. This code can output the final
         solution to a user specified output file in either double precision or
         full precision rational numbers. If you intend to use SLIP LU within
         another program, please refer to examples for help with this.

SLIPLU Followed by the listed args:

help. e.g., SLIPLU help, which indicates SLIPLU to print to guideline
for using this function.

f (or file) Filename. e.g., SLIPLU f MATRIX_NAME RHS_NAME, which indicates
SLIPLU will read matrix from MATRIX_NAME and right hand side from RHS_NAME.
For this demo, the matrix is stored in a triplet format. Please refer to
SLIP_LU/ExampleMats for some examples

c (or check). e.g., SLIPLU c, which indicates SLIPLU will check
the correctness of the solution via A*x == b.
WARNING: This is slow and should only be used for verification

p (or piv) Pivot_param. e.g., SLIPLU p 0, which inidcates SLIPLU will use
smallest pivot for pivot scheme. Other available options are listed
as follows:
       0: Smallest pivot: Default and recommended
       1: Diagonal pivoting
       2: First nonzero per column chosen as pivot
       3: Diagonal pivoting with tolerance for smallest pivot
       4: Diagonal pivoting with tolerance for largest pivot
       5: Largest pivot

q (or col) Column_order_param. e.g., SLIPLU q 0, which indicates SLIPLU
will use COLAMD for column ordering. Other available options are:
       0: None: Not recommended for sparse matrices
       1: COLAMD: Default
       2: AMD

t (or tol) tolerance_param. e.g., SLIPLU t 1e-10, which indicates SLIPLU
will use 1e-10 as the tolerance for pivot scheme 3 and 4 mentioned above.
Therefore, it is only necessary if pivot scheme 3 or 4 is used.

o2 (or out2). e.g., SLIPLU o2 1, which indicates SLIPLU will output the
errors and warnings during the process. Other available options are:
       0: print nothing
       1: just errors and warnings: Default
       2: terse, with basic stats from COLAMD/AMD and SLIP and solution

o (or out) output_param. e.g., SLIPLU o 1, which indicates SLIPLU will
output the solution as full precision rational arithmetic stdout or file
specified by the stdout redirection. It is noted that solution will be
output only when o2 >= 2. Other available options are:
       1: The solution will be output in full precision rational arithmetic
       2: The solution will be output in double precision
       3: The solution will be output in some user specified precision.
          This must be input as o 3 USER_PRECISION. Precision must be
          an integer, e.g., SLIPLU o 3 128.

If none of the above args is given, they are set to the following default:

  mat_name = "../ExampleMats/10teams_mat.txt"
  rhs_name = "../ExampleMats/10teams_v.txt"
  no solution checking
  p = 0, i.e., using smallest pivot
  q = 1, i.e., using COLAMD
  t = 0.1, not being using since p != 3 or 4
  o2 = 0, i.e., nothing will be printed
  o = 1, not being using since o2 < 2



*********example*********
Purpose: Demonstrate the use of SLIP LU for a user given input matrix

*********example2*********
Purpose: Demonstrate the use of SLIP LU for a matrix to be read in

*********example3*********
Purpose: Demonstrate the use of SLIP LU for a matrix stored in mpfr precision

*********example4*********
Purpose: Demonstrate the use of SLIP LU for multiple RHS vectors

*********example5*********
Purpose: Demonstrate the use of SLIP LU for multiple RHS vectors
         read in from file
