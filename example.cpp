#include "../Headers/SLIP_LU_config.h"

/* This example shows how to use SLIP LU with a given input matrix and a double output */
int main () {
	/* Generate a random dense 50*50 matrix */
	int n = 50, nz = 2500, num = 0;
	int* i = new int [nz];
	int* j = new int [nz];
	double* x = new double [nz];
	double** b = SLIP_initialize_double_mat(n,1);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-10,10);
	for (int k = 0; k < n; k++)
	{
		b[k][0] = distribution(generator);
		for (int p = 0; p < n; p++)
		{
			i[num] = k;
			j[num] = p;
			x[num] = distribution(generator);
			num+=1;
		}
	}
	
	/* Declare our data structures */
	SLIP_mat* A = new SLIP_mat;	        					// Input matrix
	SLIP_col* S = new SLIP_col;	        					// Column permutation 
	SLIP_LU_Options* option = new SLIP_LU_Options;		    // Factorization options
	mpz_t** b2 = SLIP_initialize_mpz_mat(n,1);          	// Right hand side vector
	double** soln = SLIP_initialize_double_mat(n,1);    	// Solution vector (double precision)
	S->q = new int [n+1];
	SLIP_Set_Options_Defaults(option);						// Set option defaults. Can change any options here
	
	/* Read in A and b */
	std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
	int type = 1;				        // Double A and b
	SLIP_create_mat_trip(n, nz, i, j, NULL, x, NULL, NULL, NULL, 0, type, A, option);	// Populate A with Triplet input
	SLIP_create_rhs(n, 1, NULL, b, NULL, NULL, NULL, 0, type, b2, option);				// Populate b
	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	option->t_inp = std::chrono::duration_cast<std::chrono::duration<float>>(t_end - t_begin);
	
	/* Factorize */
	SLIP_LU_Symbolic(A, S, option, b2, 1);
	SLIP_LU_double(A, S, b2, soln, 1, option);
	
	/* Free memory */
	delete[] i; delete[] j; delete[] x;
	SLIP_delete_double_mat(soln, n, 1);
	SLIP_delete_double_mat(b, n, 1);
	SLIP_delete_mpz_mat(b2,n,1);
	SLIP_free_memory(A, S, option);
}