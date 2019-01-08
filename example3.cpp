#include "../Headers/SLIP_LU_config.h"
/* This example uses an ccf input in mpfr precision */
int Ap[5] = {0, 3, 5, 8, 11};
int Ai[11]       = {0, 1, 2, 2, 3, 1, 2, 3, 0, 1,  2};
double Axnum[11] = {1, 2, 7, 1, 2, 4, 1, 3, 1, 12, 1};	// Numerator of x
double Axden[11] = {3, 3, 6, 1, 7, 1, 1, 1, 5, 1,  1};	// Denominator of x
double bxnum[4] = {17, 182, 61, 67};	                // Numerator of b
double bxden[4] = {15,  3,   6,  7};	                // Denominator of b

int main () {
	/* Get matrix */
	int n = 4, prec = 256, nz = 11, j;	                        // 256 bit precision
	mpfr_t ** b = SLIP_initialize_mpfr_mat(n,1,prec);
	int* i = new int [nz];
	int* p = new int [n+1];
	mpfr_t* x = SLIP_initialize_mpfr_array(nz, prec);		
	for (j = 0; j < n; j++)	                        // Get p & b
	{
		p[j] = Ap[j];
		mpfr_set_d(b[j][0], bxnum[j], MPFR_RNDN);
		mpfr_div_d(b[j][0], b[j][0], bxden[j], MPFR_RNDN);
	}
	p[n] = Ap[n];									// Finalize p
	for (j = 0; j < nz; j++)	                    // Get i and x
	{
		i[j] = Ai[j];
		mpfr_set_d(x[j], Axnum[j], MPFR_RNDN);
		mpfr_div_d(x[j], x[j], Axden[j], MPFR_RNDN);
	}
	/* Declare our data structures */
	SLIP_mat* A = new SLIP_mat;
	SLIP_col* S = new SLIP_col;
	SLIP_LU_Options* option = new SLIP_LU_Options;
	mpz_t** b2 = SLIP_initialize_mpz_mat(n,1);
	double** soln = SLIP_initialize_double_mat(n,1);
	S->q = new int [n+1];
	SLIP_Set_Options_Defaults(option);
	
	/* Read in A and b */
	std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
	int type = 4;	// Tells SLIP that we're getting MPFR input
	SLIP_create_mat_ccf(p, n, nz, i, NULL, NULL, NULL, NULL, x, prec, type, A, option);
	SLIP_create_rhs(n, 1, NULL, NULL, NULL, NULL, b, prec, type, b2, option);
	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	option->t_inp = std::chrono::duration_cast<std::chrono::duration<float>>(t_end - t_begin);
	
	/* Factorize */
	option->check = 1;	// We want to check the solution
	SLIP_LU_Symbolic(A, S, option, b2, 1);
	SLIP_LU_double(A, S, b2, soln, 1, option);

	/* Free memory */
	delete[] i; delete[] p;
	SLIP_delete_double_mat(soln, n, 1);
	SLIP_delete_mpfr_mat(b, n, 1);
	SLIP_delete_mpfr_array(x,nz);
	SLIP_delete_mpz_mat(b2,n,1);
	SLIP_free_memory(A, S, option);	
}
