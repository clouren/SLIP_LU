#ifndef SLIP_Verify
#define SLIP_Verify

/* Check the solution of the linear system. This performs a very quick rational arithmetic A*x=b
   Arguments:
   A: input matrix
   x: solution vector
   n: size of matrix
   numRHS: number of RHS vectors
   b: right hand side vector
   pinv: row permutation
   option: captures timing stats */
int SLIP_LU_Check(SLIP_mat* A, mpq_t** x, int n, int numRHS, mpz_t** b, int* pinv, SLIP_LU_Options* option)
{
	std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
	int p, j, i, check = 1;							            // Declare vars
	mpq_t temp; mpq_init(temp);
	mpq_t** b2 = SLIP_initialize_mpq_mat(n, numRHS);			// b2 stores the solution of A*x
	for (i = 0; i < numRHS; i++)
		for (j = 0; j < n; j++)
			for (p = A->p[j]; p < A->p[j+1]; p++)				// Note, there is no mpq_addmul
			{
				mpq_set_z(temp, A->x[p]);			            // temp = A[p][j]
				mpq_mul(temp, temp, x[j][i]);			        // temp = temp*x[j]
				mpq_add(b2[A->i[p]][i], b2[A->i[p]][i], temp);	// b2[p] = b2[p]-temp
			}
	for (i = 0; i < numRHS; i++)
		for (j = 0; j < n; j++)
		{
			mpq_set_z(temp, b[j][i]);				            // z = b[j] (correct b)
			if ( mpq_equal(temp, b2[j][i]) == 0)			    // Not equal
				check = 0;
		}
	mpq_clear(temp);							                // Free memory of temp and b2
	SLIP_delete_mpq_mat(b2, n, numRHS);
	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	option->t_c = std::chrono::duration_cast<std::chrono::duration<float>>(t_end - t_begin);
	return check;
}

#endif