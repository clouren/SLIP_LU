#include "../Headers/SLIP_LU_config.h"
/* This example shows how to use SLIP LU within your code and read in a matrix stored in MM format. Also shows
   how to use SLIP with an mpq output */

int main () {
	/* Declare our data structures */
	SLIP_mat* A = new SLIP_mat;
	SLIP_col* S = new SLIP_col;
	SLIP_LU_Options* option = new SLIP_LU_Options;
		
	/* Allocate memory and set defaults, read in A and b */
	SLIP_Set_Options_Defaults(option);
	option->mat_name = "../ExampleMats/10teams.mat";
	option->rhs_name = "../ExampleMats/10teams.v";
	SLIP_read_matrix(A, option->mat_name);		// Read in A
	int n = A->n;					            // Set n
	mpz_t** b = SLIP_initialize_mpz_mat(n,1);
	mpq_t** x = SLIP_initialize_mpq_mat(n,1);
	SLIP_read_rhs(b, n, option->rhs_name);
	
	/* Symbolic Ordering and Factorization */
	S->q = new int [n+1];
	SLIP_LU_Symbolic(A, S, option, b, 1);	        // Symbolic Analysis
	SLIP_LU_mpq(A, S, b, x, 1, option);        // Factorization
	
	/* Free memory */
	SLIP_free_memory(A, S, option);
	SLIP_delete_mpz_mat(b, n, 1);
	SLIP_delete_mpq_mat(x, n, 1);
}