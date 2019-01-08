#include "./Headers/SLIP_LU_config.h"

/* This program will exactly solve the sparse linear system Ax = b by performing the SLIP LU factorization.
   Please refer to README.txt for information on how to properly use this code */   
int main( int argc, char* argv[])
{
	/* Declare memory & Process Command Line */
	int n, check;
	SLIP_mat *A, *L, *U;
	SLIP_LU_Options *option = new SLIP_LU_Options;
	SLIP_Set_Options_Defaults(option);							 			 // Default options. May be changed in SLIP_LU_config.h
	if( SLIP_LU_process_command_line(argc, argv, option) == 0) return 0;	 // Process the command line

	/* Allocate memory */	
	A = new SLIP_mat; L = new SLIP_mat; U = new SLIP_mat;					// Allocate A L and U
	SLIP_read_matrix(A, option->mat_name);									// Read in the matrix specified by option->mat_name and store it in A
	n = A->n;
	mpz_t** b = SLIP_initialize_mpz_mat(n,1);
	int* pinv = new int [n];
	mpz_t* rhos = SLIP_initialize_mpz_array(n);
	SLIP_read_rhs(b, n, option->rhs_name);									// Create RHS
	
	/* Perform Column ordering */
	SLIP_col* S = new SLIP_col;
	S->q = new int [n+1];
	SLIP_LU_Symbolic(A, S, option, b, 1);									// Column ordering using either AMD, COLAMD, UMFPACK or nothing
	if (option->print == 1) SLIP_print_options(option);
	
	/* Factorization Step */
	if (SLIP_LU_Factor(A, L, U, S, rhos, pinv, option) < 0) return 0;		// SLIP LU Factorization

	/* Solve linear system */
	mpq_t** x = SLIP_initialize_mpq_mat(n,1);
	SLIP_Solve(x, b, rhos, L, U, pinv, option, 1);							// Solve LDU x = b

	/* Soln verification */
	SLIP_Permute_x(x, n, 1, S->q);											// x = Q x
	if (option->check == 1) 
		check = SLIP_LU_Check(A, x, n, 1, b, S->q, option);
	SLIP_Scale_x(x, n, 1, option);
	
	/* Output */
	SLIP_LU_Print_Stats(option, check, L->nz+U->nz-n, n, 1, x);
	
	/* Free Memory */
	SLIP_delete_mpz_mat(b,n,1); SLIP_delete_mpz_array(rhos,n); SLIP_delete_mpq_mat(x,n,1);
	SLIP_mat_delete(A); SLIP_mat_delete(L); SLIP_mat_delete(U);
	SLIP_delete_options(option); 
	SLIP_delete_col(S); 
	delete[] pinv; 
}