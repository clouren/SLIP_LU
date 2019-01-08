#ifndef SLIP_config
#define SLIP_config

/*	This software package exactly solves a sparse system of linear equations using the SLIP LU 
	factorization. 
	
	The theory associated with this software can be found in the paper (submitted to SIAM journal 
	on matrix analysis and applications):
	
	"Asymptotically Optimal Exact Solution of Sparse Linear Systems via Left-Looking 
	Roundoff-Error-Free LU Factorization"
	
	If you use this code, you must do the following:
		1) Download and install Tim Davis' Suitesparse, particularly the COLAMD and AMD routines
		   This can be obtained at http://faculty.cse.tamu.edu/davis/suitesparse.html
		2) Download and install the GMP and MPFR libraries. GMP and MPFR can be found at
		   https://gmplib.org/
		   http://www.mpfr.org/
		   
	If you use it for a publication, please cite the associated paper
	
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Authors--------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

		Christopher Lourenco, Erick Moreno-Centeno, Adolfo Escobedo, and Timothy Davis

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Copyright------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

	This software is licensed under the GNU General Public License
	See license.txt for license info.

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------DISCLAIMER-----------------------------------------------------------------
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

********** This software comes with no implied warranty, use it at your own risk *********************

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Summary--------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
	
	This software package solves the linear system Ax = b exactly. The input matrix and right
	hand side vectors are stored as either integers, double precision numbers, multiple precision
	floating points (through the mpfr library) or as rational numbers (as a collection of numerators
	and denominators using the GMP mpq_t data structure). Appropriate routines within the code transform
	the input into an integral matrix in compressed column form.
	
	This package computes the factorization PAQ = LDU. Note that we store the "functional" 
	form of the factorization by only storing L and U. The user is given some freedom to select
	the permutation matrices P and Q. The recommended default settings select Q using the COLAMD column
	ordering and selects P via a partial pivoting scheme in which the smallest entry in 
	column k is selected as the kth pivot. Alternative strategies allowed to select Q include the AMD
	column ordering, no column permutation (Q=I), or using the UMFPACK P and Q. For pivots, there are a 
	variety of potential schemes including traditional partial pivoting, diagonal pivoting, tolerance pivoting etc. 
	This package does not allow pivoting based on sparsity criterion.
	
	The factors L and U are computed via integer preserving operations via integer-preserving
	Gaussian elimination. The key part of this algorithm is a REF Sparse triangular solve 
	function which exploits sparsity to reduce the number of operations that must be performed.
	
	Once L and U are computed, a simplified version of the triangular solve is performed which assumes
	the vector b is dense. The final solution vector x is gauranteed to be exact. This vector can be 
	output in one of three ways: 1) full precision rational arithmetic (as a sequence of numerators 
	and denominators) using the GMP mpq_t data type, 2) double precision while not exact will produce
	a solution accurate to machine roundoff unless the size of the associated solution exceeds double
	precision (i.e., the solution is 10^500 or something), 3) variable precision floating point using 
	the GMP mpfr_t data type. The associated precision is user defined.
	

------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------C++ Libraries--------------------------------------------------------------
------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------*/

/* Standard C/C++ libraries */

# include <stdlib.h>
# include <iostream>
# include <fstream>
# include <chrono>
# include <math.h>
# include <sstream>
# include <random>
# include <iomanip>
# include <string>
# include <algorithm>
# include <new> 

/* GMP and MPFR headers */
# include <gmp.h>
# include <gmpxx.h>
# include <mpfr.h>

/* SuiteSparse headers */
# include "SuiteSparse_config.h"
# include "colamd.h"    
# include "amd.h"   

/*----------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Default Parameters---------------------------------------------------------
------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------*/

/* Current version of the code */
#define SLIP_LU_Version "0.4"

/* Primary author of code */
#define Author "Christopher Lourenco, Erick Moreno-Centeno, Adolfo Escobedo, Timothy Davis"

/* Tolerance used in the pivoting schemes. This number can be anything in between 0 and 1.
   A value of 0 selects the diagonal element exclusively and a value of 1 selects the smallest or
   largest pivot exclusively only in a tolerance pivoting based method */
#define DEFAULT_TOL 0.1

/* Default matrix and right hand side used for 'SLIP_LU.cpp' */
#define DEFAULT_MAT "./ExampleMats/10teams.mat"
#define DEFAULT_RHS "./ExampleMats/10teams.v"

/* Check parameter. If this = 1 then the solution to the system is checked for accuracy */
#define DEFAULT_CHECK 0

/* Pivoting scheme used for SLIP LU. 0: Smallest, 1: Natural/Diagonal pivoting 2: Choose first nonzero
   3: Diagonal with tolerance and smallest pivot 4: Diagonal with tolerance and largest pivoting
   5: Largest pivot */
#define DEFAULT_PIVOT 0

/* Column ordering used. 0: COLAMD, 1: AMD, 2: no column ordering */
#define DEFAULT_ORDER 0

/* Defines printing to be done */
#define DEFAULT_PRINT 0
#define DEFAULT_PRINT2 0
#define DEFAULT_PRINT3 1
#define DEFAULT_OUTFILE "default_output.out"
#define DEFAULT_RAT 1

/* MPFR precision used (quad is default) */
#define DEFAULT_PRECISION 128

/* Needed for MATLAB mex files */
#ifdef MATLAB
	# define INT ptrdiff_t
	# include "mex.h"
#endif


/* For random matrices to compare to REF LU */

#define DEFAULT_N 100
#define DEFAULT_NUMRHS 50
#define DEFAULT_LB -99
#define DEFAULT_UB 99
#define DEFAULT_SEED1 11
#define DEFAULT_SEED2 1400
#define DEFAULT_DENSITY 1.0

/*----------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Common Macros--------------------------------------------------------------
------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------*/

/* Used for the reach functions. Taken from CSparse */
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }

/*----------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Data Structures------------------------------------------------------------
------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------*/

/* This struct defines the command line options and some timing statistics for the factorization. In
   addition, it stores the determinant of the matrix and scaling parameter used on A and b. In essence
   it serves as a global struct to define all options */
typedef struct SLIP_LU_Options	
{
	int check;				            // 1 If the solution to the linear system will be checked
	int pivot;			                // Type of pivoting scheme used.
	int order;			               	// 0: COLAMD, 1: AMD, 2: None, 3: UMFPACK >3 or < 0 defaults to COLAMD
	int print;          				// 0: Nothing printed, 1: COLAMD/AMD and SLIP stats printed 
	int print2;			            	// 0: output file not used, 1: output file used
	int print3;         				// 0: timing stats not printed, 1: timing stats printed
	int rat;		              		// 1 if output file is to print in rational, 2 for double, 3 for MPFR
	int prec;			             	// Precision used to output file if MPFR is chosen
	double tol;				            // User specified tolerance
	std::string mat_name;			    // Name of matrix stored in matrix market format
	std::string rhs_name;			    // Name of RHS vector stored as a single column vector
	std::string out_file;			    // Output file if desired
	std::chrono::duration<float> t_sym;	// Symbolic analysis time
	std::chrono::duration<float> t_inp;	// Input manipulation time
	std::chrono::duration<float> t_c;	// Solution check time (if necessary)
	std::chrono::duration<float> t_f;	// LU Factorization time
	std::chrono::duration<float> t_s;	// Forward/back solve time
	std::chrono::duration<float> t_tot;	// Total run time	
	mpq_t LU_scale;				        // Parameter used to scale the A matrix
	mpq_t b_scale;			        	// Parameter used to scale the right hand side vector
	mpz_t determinant;	        		// Determinant of the matrix
} SLIP_LU_Options;

/* This struct defines a matrix stored in sparse compressed column form. Since this code deals with exact data structures,
   each SLIP_mat stores both the size of internal vectors and the number of elements allocated in each one. This is done
   to use as little memory as possible */
typedef struct SLIP_mat
{
	int n;		// Number of columns
	int m;		// Number of rows. Code assumes n = m
	int nzmax;	// Allocated size of A->i and A->x
	int nz;		// Number of nonzeros in the matrix 
	int* p;		// Column pointers. Array size is n+1
	int* i;		// Row indices. Array size is nzmax, number of entries is nz
	mpz_t* x;	// Values in matrix. Array size is nzmax, number of entries is nz
} SLIP_mat;

/* This struct stores the column permutation for LU and the guess on nnz for L and U */
typedef struct SLIP_col
{
	int* q;		// Fill reduced column permutation for LU
	int lnz;	// Approximate number of nonzeros in L. i.e., initial size of L matrix
	int unz;	// Approximate number of nonzeros in U. i.e., initial size of U matrix
} SLIP_col;

/* This struct stores a sparse integral matrix in triplet form. This structure is only
   used if the matrix is read in from a file in integral matrix market format */
typedef struct SLIP_trip		
{
	int n;		// Number of rows
	int m;		// Number of columns
	int* i;		// Row pointers. Size nz
	int* j;		// Column pointers. Size nz
	mpz_t* x;	// Values in matrix. Size nz
	int nz;		// Number of nonzeros in matrix
} SLIP_trip;

/*----------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------
---------------------------Primary Routines-----------------------------------------------------------
------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------*/

/* Purpose: This code prints information about SLIP LU.
   Arguments: None  */
void SLIP_LU_Info()
{
	std::cout<<"\n****Software Information****";
	std::cout<<"\n\nYou are using SLIP_LU Version: "<<SLIP_LU_Version;
	std::cout<<"\nThis code accompanies the paper:\n\n\t"<<PAPER<<"\n";
	std::cout<<"\nThis code is copyright by "<<Author <<"\n\n";
}

/* Purpose: This function initializes an mpz array of size n and allocates default size.
   Arguments:
   n: Size of x */
mpz_t* SLIP_initialize_mpz_array (int n)		
{
	mpz_t* x = new mpz_t [n];
	for (int i = 0; i <n; i++)
		mpz_init(x[i]);
	return x;
}

/* Purpose: This function initializes an mpz array of size n and allocates memory for
   numbers of bit size prec. If the relative size of numbers is known ahead of time, this 
   is more efficient than the SLIP_initialize_mpz_array
   Arguments:
   n: size of the array
   size: Relative size of numbers */
mpz_t* SLIP_initialize_mpz_array2(int n, int size)
{
	mpz_t* x = new mpz_t [n];
	for (int i = 0; i < n; i++)
		mpz_init2(x[i], size);
	return x;
}

/* Purpose: This function initializes an mpq array of size n. 
   This function must be called for all mpq arrays created.
   Arguments:
   n: size of the array */
mpq_t* SLIP_initialize_mpq_array (int n)		
{
	mpq_t* x = new mpq_t [n];
	for (int i = 0; i < n; i++)
		mpq_init(x[i]);
	return x;
}

/* Purpose: This function initializes an double array of size n and sets each element equal to zero
   Arguments:
   n: size of x */
double* SLIP_initialize_double_array (int n)
{
	double* x = new double [n];
	for (int i = 0; i < n; i++)
		x[i] = 0;
	return x;
}

/* Purpose: This function initializes an integer array of size n and sets each element equal to zero
   Arguments:
   n: size of x */
int* SLIP_initialize_int_array (int n)
{
	int* x = new int [n];
	for (int i = 0; i < n; i++)
		x[i] = 0;
	return x;
}

/* Purpose: This function initializes a MPFR array of desired precision
   Arguments:
   n: size of the array
   prec: associated floating point precision */
mpfr_t* SLIP_initialize_mpfr_array( int n, int prec)
{
	mpfr_t* x = new mpfr_t [n];
	for (int i = 0; i < n; i++)
		mpfr_init2(x[i], prec);
	return x;
}

/* Purpose: This function initializes an int matrix of size m*n
   Arguments:
   m: number of rows
   n: number of columns */
int** SLIP_initialize_int_mat(int m, int n)
{
    int **x;
    x = new int *[m];
    for (int j = 0; j < m; j++)
        x[j] = new int [n];
    return x;
}

/* Purpose: This function initializes a double matrix of size m*n
   Arguments:
   m: number of rows
   n: number of columns */
double** SLIP_initialize_double_mat(int m, int n)
{
	double **x;
	x = new double *[m];
	for (int j = 0; j < m; j++)
		x[j] = new double [n];
	return x;
}

/* Purpose: This function initializes a dense mpz_t matrix of size m*n to default size
   Arguments:
   m: number of rows
   n: number of columns */
mpz_t** SLIP_initialize_mpz_mat(int m, int n)
{
	mpz_t **x;
	x = new mpz_t *[m];
	for (int j = 0; j < m; j++)
	{
		x[j] = new mpz_t [n];
		for (int i = 0; i < n; i++)
			mpz_init(x[j][i]);
	}
	return x;
}

/* Purpose: This function initializes a mpq_t matrix of size m*n
   Arguments:
   m: number of rows
   n: number of columns */
mpq_t** SLIP_initialize_mpq_mat(int m, int n)
{
	mpq_t **x;
	x = new mpq_t *[m];
	for (int j = 0; j < m; j++)
	{
		x[j] = new mpq_t [n];
		for (int i = 0; i < n; i++)
			mpq_init(x[j][i]);
	}
	return x;
}

/* Purpose: This function initializes a mpfr_t matrix of size m*n with precision prec
   Arguments:
   m: number of rows
   n: number of columns
   prec: floating point precision */
mpfr_t** SLIP_initialize_mpfr_mat(int m, int n, int prec)
{
	mpfr_t **x;
	x = new mpfr_t *[m];
	for (int j = 0; j < m; j++)
	{
		x[j] = new mpfr_t [n];
		for (int i = 0; i < n; i++)
			mpfr_init2(x[j][i],prec);
	}
	return x;
}

/* Purpose: This function resets an mpz array of size n by setting each term equal to zero.
   This is more efficient then deletion and reallocation/reinitialization
   Arguments:
   x: mpz array to be reset
   n: size of x */
void SLIP_reset_mpz_array(mpz_t* x, int n)
{
	for (int i = 0; i < n; i++)
		mpz_set_ui(x[i],0);
}

/* Purpose: This function resets an mpz array of size n with the nonzero pattern given.
   This is more efficient than iterating accross all nonzeros in vector x. 
   Arguments:
   x: mpz array to be reset
   n: size of x
   top: beginning of the nonzero pattern
   xi: nonzero pattern */
void SLIP_reset_mpz_array_2(mpz_t* x , int n, int top, int* xi)
{
	for (int i = top; i < n; i++)	// Access the nonzero pattern located in xi[top..n-1]
		mpz_set_ui(x[xi[i]], 0);	// Each nonzero x[i] = 0
}

/* Purpose: This function initializes an int vector of size n and sets the value equal to -1.
   This function is used for the history and pivot vectors. 
   Arguments: 
   h: int vector to be reset
   n: size of the int vector */
void SLIP_reset_int_array(int* h, int n)	
{
	for (int i = 0; i < n; i++)
		h[i] = -1;
}

/* Purpose: This function resets an int vector of size n and sets each term equal to -1 with nonzero pattern given.
   This is more efficient than resetting each term individually
   Arguments:
   h: int vector to be reset
   n: size of h
   top: beginning of nonzero pattern
   xi: nonzero pattern */
void SLIP_reset_int_array_2 (int* h, int n, int top, int* xi)
{
	for (int i = top; i < n; i++)	// Access the nonzero pattern located in xin[top..n-1]
		h[xi[i]] = -1;		        // Each "nonzero" h[i] = -1
}

/* Purpose: This function clears the memory used for an mpz vector of size n. 
   Call this function for all mpz vectors when done.
   Arguments:
   x: mpz array to be deleted
   n: Size of x */
void SLIP_delete_mpz_array( mpz_t* x, int n)
{
	for (int i = 0; i < n; i++)
		mpz_clear(x[i]);
	delete[] x;
}

/* Purpose: This function clears the memory used for an mpq vector of size n. 
   Call this for all mpq vectors when done.
   Arguments:
   x: mpq array to be deleted
   n: size of x */
void SLIP_delete_mpq_array( mpq_t* x, int n)
{
	for (int i = 0; i < n; i++)
		mpq_clear(x[i]);
	delete[] x;
}

/* Purpose: This function clears the memory used for an mpfr array of size n
   Arguments:
   x: mpfr array to be deleted 
   n: size of x */
void SLIP_delete_mpfr_array(mpfr_t* x, int n)
{
	for (int i = 0; i < n; i++)
		mpfr_clear(x[i]);
	delete[] x;
}

/* Purpose: This function deletes a dense mpfr matrix
   Arguments:
   A: Dense mpfr matrix
   n: number of columns of A
   m: number of rows of A */
void SLIP_delete_mpfr_mat(mpfr_t **A, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			mpfr_clear(A[i][j]);
		delete[] A[i];
	}
	delete[] A;
}

/* Purpose: This function deletes a dense mpz matrix
   Arguments:
   A: The dense mpz matrix
   n: number of columns of A
   m: number of rows of A */
void SLIP_delete_mpz_mat(mpz_t**A, int m, int n)
{
	for (int i = 0; i < m; i++)			// Iterate accross all entries and clear the individual memory
	{
		for (int j = 0; j < n; j++)
			mpz_clear(A[i][j]);
		delete[] A[i];
	}
	delete[] A;
}

/* Purpose: This function deletes a dense mpq matrix
   Arguments:
   A: dense mpq matrix
   n: number of columns of A
   m: number of rows of A */
void SLIP_delete_mpq_mat(mpq_t**A, int m, int n)
{
	for (int i = 0; i < m; i++)			// Iterate accross all entries and clear the individual memory
	{
		for (int j = 0; j < n; j++)
			mpq_clear(A[i][j]);
		delete[] A[i];
	}
	delete[] A;
}

/* Purpose: This function deletes a dense double matrix.
   Arguments:
   x: dense matrix
   n: number of columns
   m: number of rows */
void SLIP_delete_double_mat(double** A, int m, int n)
{
	for (int i = 0; i < m; i++)
		delete[] A[i];
	delete[] A;
}

/* Purpose: This function deletes a dense int matrix.
   Arguments:
   x: dense matrix
   n: number of columns
   m: number of rows */
void SLIP_delete_int_mat(int** A, int m, int n)
{
	for (int i = 0; i < m; i++)
		delete[] A[i];
	delete[] A;
}

/* Purpose: This function allocates a SLIP LU matrix of size n*m with array size nzmax
   This function manually initializes each entry in A->x therefore they are immediately
   ready to be operated on. This is less efficient but more user friendly.
   Arguments:
   n: number of columns
   m: number of rows (recall m=n assumed)
   nzmax: size of allocated i and x arrays
   A: sparse matrix data structure to be allocated */
void SLIP_mat_alloc (int n, int m, int nzmax, SLIP_mat* A)
{
	A->m = m;			                		// Rows of the matrix
	A->n = n;			                		// Columns of the matrix
	A->nz = 0;				                	// Currently 0 nonzeros
	A->nzmax = nzmax;				            // Size of the vectors
	A->x = SLIP_initialize_mpz_array(nzmax);	// Initialize the mpz array x
	A->p = SLIP_initialize_int_array(n+1); 		// Initialize p
	A->i = SLIP_initialize_int_array(nzmax); 	// Initialize i
}

/* Purpose: This function allocates a SLIP LU matrix of size n*m with array size nzmax
   This version does not allocate individual the values in x. As a result, it is more
   memory efficient, but also less user friendly.
   Arguments:
   n: number of columns
   m: number of rows (recall m=n assumed)
   nzmax: size of allocated i and x arrays
   A: sparse matrix data structure to be allocated */
void SLIP_mat_alloc2 (int n, int m, int nzmax, SLIP_mat* A)
{   
	A->m = m;					                // Rows of the matrix
	A->n = n;					                // Columns of the matrix
	A->nz = 0;					                // Currently 0 nonzeros
	A->nzmax = nzmax;				            // Size of the vectors
	A->x = new mpz_t [nzmax];			        // Allocate memory for x values
	A->p = SLIP_initialize_int_array(n+1); 		// Initialize p
	A->i = SLIP_initialize_int_array(nzmax); 	// Initialize i
}

/* Purpose: This function expands a SLIP LU matrix by doubling its size
   Arguments:
   A: the matrix to be expanded */
void SLIP_mat_realloc ( SLIP_mat* A)
{
	mpz_t* temp_x = SLIP_initialize_mpz_array(2*A->nzmax);	// Initialize temp x
	int *  temp_i = SLIP_initialize_int_array(2*A->nzmax);	// Initialize temp r
	
	for (int i = 0; i < A->nz; i++)			        // Copy over the old values
	{
		temp_i[i] = A->i[i];
		mpz_swap(temp_x[i],A->x[i]);		        // temp_x[i] = A->x[i], A->x[i] = temp_x[i]
	}
	SLIP_delete_mpz_array(A->x,A->nzmax);		    // Delete old x
	delete[] A->i;					                // Delete old i
	A->x = temp_x;					                // Set new x
	A->i = temp_i;					                // Set new r
	A->nzmax = A->nzmax*2; 
}

/* Purpose: This function expands a SLIP LU matrix by doubling its size. This version merely expands
   x and i and does not initialize/allocate the values!
   Arguments:
   A: the matrix to be expanded */
void SLIP_mat_realloc2 ( SLIP_mat* A)
{
	int size;				                				// Temporary pointers that will be the new vectors
	mpz_t* temp_x = new mpz_t [2*A->nzmax];		        	// Allocate memory for temp x
	int*   temp_i = SLIP_initialize_int_array(2*A->nzmax);	// Initialize temp i
	for (int i = 0; i < A->nz; i++)		        			// Copy over the old values
	{
		temp_i[i] = A->i[i];
		size = mpz_sizeinbase(A->x[i], 2);
		mpz_init2(temp_x[i], size+2);
		mpz_swap(temp_x[i],A->x[i]);		        // temp_x[i] = A->x[i], A->x[i] = temp_x[i]
	}
	SLIP_delete_mpz_array(A->x,A->nz);		        // Delete old x
	delete[] A->i;					                // Delete old i
	A->x = temp_x;					                // Set new x
	A->i = temp_i;					                // Set new i
	A->nzmax = A->nzmax*2; 
}

/* Purpose: This function collapses a SLIP matrix. Essentially it shrinks the size of x and i 
   so that they only take up the number of elements in the matrix. 
   For example if A->nzmax = 1000 but A->nz = 500, r and x are of size 1000 
   so this function shrinks them to size 500
   Arguments:
   A: matrix to be shrunk */
void SLIP_mat_collapse (SLIP_mat* A)
{
	int*   temp_i = new int [A->nz];			        // Declare size of r
	mpz_t* temp_x = SLIP_initialize_mpz_array(A->nz);	// Initialize x
	
	for (int i = 0; i < A->nz; i++)			    // Copy over values
	{
		temp_i[i] = A->i[i];
		mpz_set(temp_x[i],A->x[i]);
	}
	SLIP_delete_mpz_array(A->x,A->nzmax);		// Delete memory associated with x and r
	delete[] A->i;
	A->x = temp_x;
	A->i = temp_i;					            // Assign new values
	A->nzmax = A->nz;
}

/* Purpose: This function collapses a SLIP matrix. Essentially it shrinks the size of x and i 
   so that they only take up the number of elements in the matrix. 
   For example if A->nzmax = 1000 but A->nz = 500, r and x are of size 1000 
   so this function shrinks them to size 500
   NOTE: This differs from the above function by the assumption of the size of A->x. This version
   assumes that we are performing the "memory efficient" version of the code.
   Arguments:
   A: matrix to be shrunk */
void SLIP_mat_collapse2 (SLIP_mat* A)
{
	int*   temp_i = new int [A->nz];			        // Declare size of r
	mpz_t* temp_x = SLIP_initialize_mpz_array(A->nz);	// Initialize x
	
	for (int i = 0; i < A->nz; i++)			    // Copy over values
	{
		temp_i[i] = A->i[i];
		mpz_set(temp_x[i],A->x[i]);
	}
	SLIP_delete_mpz_array(A->x,A->nz);		    // Delete memory associated with x and r
	delete[] A->i;
	A->x = temp_x;
	A->i = temp_i;					            // Assign new values
	A->nzmax = A->nz;
}


/* Purpose: This function deletes the sparse matrix A
   Arguments: 
   A: matrix to be deleted */
void SLIP_mat_delete(SLIP_mat* A)
{
	SLIP_delete_mpz_array(A->x,A->nzmax);
	delete[] A->i; delete[] A->p; delete A;
}

/* Purpose: This function deletes the option struct
   Arguments:
   option: Struct to be deleted */
void SLIP_delete_options(SLIP_LU_Options* option)
{
	mpq_clear(option->LU_scale); mpq_clear(option->b_scale); mpz_clear(option->determinant);
	delete option;
}

/* Purpose: This function deletes the SLIP_col struct
   Arguments:
   S: Structure to be deleted */
void SLIP_delete_col(SLIP_col* S)
{
	delete[] S->q;
	delete S;
}

/* Purpose: This function frees memory of the 3 main structs used 
   preventing the user from having to call each one seperately
   Arguments:
   A: input matrix
   S: SLIP_col structure
   option: options struct */
void SLIP_free_memory(SLIP_mat* A, SLIP_col* S, SLIP_LU_Options* option)
{
	SLIP_delete_col(S);
	SLIP_mat_delete(A);
	SLIP_delete_options(option);
}

/* Purpose: This function prints onto the screen a dense mpz matrix. 
   Primarily used for debugging/troubleshooting
   Arguments:
   A: dense mpz matrix to be printed 
   n: number of columns in A
   m: number of rows in A */
void SLIP_print_mpz_mat (mpz_t** A, int n, int m)
{
	for (int i = 0; i <m; i++)
	{
		for (int j = 0; j < n; j++)
			std::cout << A[i][j] << " ";
		std::cout << "\n";
	}
}

/* Purpose: This function converts a matrix from ccf to dense with no
   permutation applied. Primarily used for debugging/troubleshooting
   Arguments:
   A: sparse matrix to be converted
   B: dense matrix to hold A */
void SLIP_convert_mpz_dense(SLIP_mat* A, mpz_t** B)
{
	int counter = 0;
	for (int i = 0; i < A->n; i++)
		for (int j = 0; j <  A->p[i+1] - A->p[i]; j++)
		{
			mpz_set(B[A->i[counter]][i],A->x[counter]);	// Row index of entry in column i
			counter+=1;
		}
}

/* Purpose: This function prints onto the screen an ccf matrix with no row or column permutation. 
   Primarily used for debugging 
   Arguments:
   A: matrix to be printed*/
void print_sp_mat (SLIP_mat* A)
{	
	int n, m;
	n = A->n;			            // Number of columns
	m = A->m;			            // Number of rows	
	mpz_t** B = SLIP_initialize_mpz_mat(m,n);
	SLIP_convert_mpz_dense(A, B);	// Get A into dense
	SLIP_print_mpz_mat(B,n,m);	    // Print out matrix
	SLIP_delete_mpz_mat(B,n,m);	    // Delete memory occupied by dense matrix
} 

#include "SLIP_LU_Factor.h"
#include "SLIP_Conversions.h"
#include "SLIP_Soln_Verify.h"
#include "SLIP_LU_input.h"
#include "SLIP_LU_Output.h"
#include "SLIP_LU.h"

#endif
