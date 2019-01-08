#include "../Headers/SLIP_LU_config.h"
/* This example shows how to use multiple RHS vectors */
int Ap[5] = {0, 3, 5, 8, 11};
int Ai[11]       = {0, 1, 2, 2, 3, 1, 2, 3, 0, 1,  2};
double Axnum[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};	// Numerator of x
double Axden[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};	// Denominator of x
double bxnum[4] = {17, 182, 61, 67};	                // Numerator of b
double bxden[4] = {15,  3,   6,  7};	                // Denominator of b
double bxnum2[4] = {8, 50, 25, 23};	                    // Numerator of b2
double bxden2[4] = {15, 3, 6, 7};	                    // Denominator of b2

int main () {
	/* Get matrix */
	int n = 4, nz = 11, numRHS = 2;
	double* x = new double [nz];
	double **b = SLIP_initialize_double_mat(n, numRHS);
	int* i = new int [nz];
	int* p = new int [n+1];
	for (int j = 0; j < n; j++)	        // Get b & p
	{
		p[j] = Ap[j];
		b[j][0] = bxnum[j]/bxden[j];
		b[j][1] = bxnum2[j]/bxden2[j];
	}
	p[n] = Ap[n];
	for (int j = 0; j < nz; j++)	    // Get i and x
	{
		i[j] = Ai[j];
		x[j] = Axnum[j]/Axden[j];
	}
	/* Declare our data structures */
	SLIP_mat* A = new SLIP_mat;
	SLIP_col* S = new SLIP_col;
	SLIP_LU_Options* option = new SLIP_LU_Options;
	mpz_t** b2 = SLIP_initialize_mpz_mat(n, numRHS);;
	double** soln = SLIP_initialize_double_mat(n, numRHS);;
	S->q = new int [n+1];
	SLIP_Set_Options_Defaults(option);
	
	/* Read in A and b */
	std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
	int type = 1;	// Tells SLIP it's getting double input
	SLIP_create_mat_ccf(p, n, nz, i, NULL, x, NULL, NULL, NULL, 0, type, A, option);
	SLIP_create_rhs(n, numRHS, NULL, b, NULL, NULL, NULL, 0, type, b2, option); 
	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	option->t_inp = std::chrono::duration_cast<std::chrono::duration<float>>(t_end - t_begin);
		
	/* Factorize */
	option->check = 1;
	SLIP_LU_Symbolic(A, S, option, b2, numRHS);
	SLIP_LU_double(A, S, b2, soln, numRHS, option);
		
	/* Free memory  */
	delete[] i; delete[] p; delete[] x;
	SLIP_delete_double_mat(b, n, numRHS);
	SLIP_delete_double_mat(soln, n, numRHS);
	SLIP_delete_mpz_mat(b2,n, numRHS);
	SLIP_free_memory(A, S, option);	
}
