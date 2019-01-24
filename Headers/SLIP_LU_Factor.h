#ifndef SLIP_Factor
#define SLIP_Factor

/* This header contains the routines and functions used for LU factorization and forward/back solve

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------LU Factorization Routines-------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* Purpose: This function obtains column k from matrix A and stores it in the dense vector x
   Arguments:
   A: input matrix
   k: column to extract
   x: A(:,k) */
void SLIP_get_column(SLIP_mat* A, int k, mpz_t* x)
{
	for (int i = A->p[k]; i < A->p[k+1]; i++) 	// Iterating accross the nonzeros in column k
		mpz_set(x[A->i[i]],A->x[i]);		    // Value of the ith nonzero
}

/* Purpose: This function performs a depth first search of the graph of the matrix starting at node j
   The output of this function is the set of nonzero indices 
   Arguments:
   j: What node to start DFS at
   L: matrix which represents the Graph of L 
   top: beginning of stack
   xi: the nonzero pattern
   pstack: workspace vector
   pinv: row permutation */
int cs_dfs (int j, SLIP_mat* L, int top, int* xi, int* pstack, int* pinv)
{
	int i, p, p2, done, jnew, head = 0;
	xi[0] = j;					        // Initialize the recursion stack
	while (head >= 0)
	{
		j = xi[head];				    // The j value of the nonzero
		jnew = pinv[j];				    // The relative j value 
		if (!CS_MARKED (L->p,j))		// Mark L->p[j] if not marked yet
		{
			CS_MARK(L->p,j);
			pstack[head] = (jnew < 0) ? 0 : CS_UNFLIP(L->p[jnew]);
		}
		done = 1;				        // Node j is done if no unvisited neighbors
		p2 = (jnew < 0) ? 0 : CS_UNFLIP(L->p[jnew+1]);
		for (p = pstack[head]; p < p2; p++) 	// Examine all neighbors of j
		{
			i = L->i[p];			    // Looking at neighbor node i
			if (CS_MARKED(L->p,i))		// Skip already visited node
				continue;
			pstack[head] = p;		    // pause DFS of node j
			xi[++head] = i;			    // Start DFS at node i
			done = 0;			        // node j is not done
			break;				        // break to start dfs i
		}
		if (done != 0)
		{
			head--;
			xi[--top] = j;
		}
	}
	return top;
}

/* Purpose: This function computes the reach of column k of A on the graph of L
   mathematically that is: xi = Reach(A(:,k))_G_L
   Arguments:
   L: matrix representing graph of L
   A: input matrix
   k: column of A of interest
   xi: nonzero pattern
   pinv: row permutation */
int cs_reach (SLIP_mat* L, SLIP_mat* A, int k, int* xi, int* pinv)
{
	int p, n = L->n, top = n;
	for (p = A->p[k]; p < A->p[k+1]; p++) 	// Iterating across number of nonzero in column k
		if (!CS_MARKED(L->p, A->i[p])) 	    // DFS at unmarked node i
			top = cs_dfs(A->i[p], L, top, xi, xi + n, pinv);
	
	for ( p = top; p < n; p++)  		    // Restore L
		CS_MARK(L->p, xi[p]);
	return top;
}

/* Purpose: This function selects the pivot element as the largest in the column
   This is activated if the user sets option->pivot = 5
   This pivoting scheme is NOT recommended for SLIP LU
   Arguments:
   x: kth column of L and U
   pivs: vector which indicates whether each row has been pivotal
   n: dimension of problem
   top: nonzero pattern is located in xi[top..n-1] 
   xi: nonzero pattern of x */
int SLIP_get_largest_pivot (mpz_t* x, int* pivs, int n, int top, int* xi )
{
	int i, inew, piv = -1;
	mpz_t big; mpz_init(big);
	for (i = top; i < n; i++)					            // Iterate accross the nonzeros in x
	{
		inew = xi[i];						                // Location of the ith nonzero
		if (pivs[inew] < 0 && mpz_cmpabs(big,x[inew]) < 0)	// inew can be pivotal
		{
			piv = inew;					                    // Current largest pivot location
			mpz_set(big,x[inew]);				            // Current largest pivot value
		}
	}
	mpz_clear(big);							                // Frees the memory occupied by the pivot value
	return piv;							                    // Return the pivot
}

/* Purpose: This function selects the pivot element that is the smallest in the column 
   This is activated by default or if the user sets option->pivot = 0
   This is the recommended pivoting scheme for SLIP LU
   Arguments:
   x: kth column of L and U
   pivs: vector indicating whether each row has been pivotal
   n: dimension of problem
   top: nonzeros are stored in xi[top..n-1]
   xi: nonzero pattern of x */
int SLIP_get_smallest_pivot (mpz_t* x, int* pivs, int n, int top, int* xi)
{
	int i, inew, j = n, flag = top, piv = -1;			// Flag is non-negative until we have an initial starting value for small
	mpz_t small; mpz_init(small);
	while (flag > -1 && flag < n)					    // Find an initial pivot. Fails if all terms are 0 in array x
	{
		inew = xi[flag];					            // i location of first nonzero
		if (pivs[inew] < 0 && mpz_sgn(x[inew]) != 0)	// inew can be pivotal
		{
			mpz_set(small,x[inew]);				        // Current smallest pivot
			piv = inew;					                // Current smallest pivot location
			j = flag;					                // Where to start the search for rest of nonzeros
			flag = -5;					                // Exit the while loop
		}
		flag += 1;						                // Increment to next nonzero to search
	}
	for (i = j; i < n; i++)						        // Iterate across remaining nonzeros
	{
		inew = xi[i];
		if (pivs[inew] < 0 && mpz_cmpabs(small,x[inew]) > 0)	// inew can be pivotal
		{
			if ( mpz_sgn(x[inew]) != 0)
			{
				piv = inew;				                // Current best pivot location
				mpz_set(small,x[inew]);			        // Current best pivot value
			}
		}
	}
	mpz_clear(small);						            // Clear memory for small
	return piv;
}

/* This function obtains the first eligible nonzero pivot
   This is enabled if the user sets option->pivot = 2
   This pivoting scheme is not recommended 
   Arguments:
   x: kth column of L and U
   pivs: vector indicating which rows are pivotal
   n: size of x
   top: nonzero pattern is located in xi[top..n-1]
   xi: nonzero pattern of x */
int SLIP_get_nonzero_pivot(mpz_t* x, int* pivs, int n, int top, int* xi)
{
	int i, inew, piv = -1;
	for (i = top; i < n; i++)				            // Iterate accross the nonzeros in x
	{
		inew = xi[i];					                // inew is the location of the ith nonzero
		if (mpz_sgn(x[inew]) != 0 && pivs [inew] < 0)	// x[inew] is an eligible pivot
		{
			piv = inew;
			break;					                    // End the loop
		}
	}
	return piv;
}

/* This function performs the pivoting for the SLIP LU factorization. The options are:
   Order is:
	0: Smallest pivot (default)
	1: Natural/Diagonal pivoting
	2: Choose first nonzero
	3: Diagonal with tolerance and smallest pivot
	4: Diagonal with tolerance and largest pivoting
	5: Largest pivot
   Arguments:
   x: kth column of L and U
   pivs: vecor indicating which rows have been pivotal
   n: dimension
   top: nonzero pattern is located in xi[top..n-1]
   xi: nonzero pattern of x
   order: tells what kind of pivoting to use 
   col: what column of A we are currently in (real kth column i.e., q[k])
   k: iteration of the algorithm
   rhos: vector of pivots
   pinv: row permutation
   row_perm: opposite of pinv. if pinv[i] = j then row_perm[j] = i 
   tolerance: tolerance used if some tolerance based pivoting is used */
int SLIP_get_pivot(mpz_t* x, int* pivs, int n, int top, int* xi, int order, int col, int k, mpz_t* rhos, int* pinv, int* row_perm, double tolerance)
{
	int pivot;
	if (order == 0)								            // Smallest pivot
		pivot = SLIP_get_smallest_pivot(x, pivs, n, top, xi);
	else if (order == 1) 							        // Diagonal
	{
			
		if ( mpz_sgn(x[col]) != 0 && pivs[col] < 0)			// Check if x[col] is eligible
			pivot = col;
		else								                // If not eligible, take smallest pivot
			pivot = SLIP_get_smallest_pivot(x, pivs, n, top, xi);
	}
	else if (order == 2)							        // First nonzero
		pivot = SLIP_get_nonzero_pivot(x, pivs, n, top, xi);
	else if (order == 3)							        // Tolerance with smallest pivot
	{
		pivot = SLIP_get_smallest_pivot(x, pivs, n, top, xi);
		if (mpz_sgn(x[col]) != 0 && pivs[col] < 0)			// Checking x[col] vs smallest pivot
		{
			mpq_t tol, ratio;
			mpq_inits(tol,ratio,NULL);				        // Initialize tolerance and ratio
			mpq_set_d(tol,tolerance);				        // Set user specified tolerance
			mpq_set_num(ratio,x[pivot]);				    // ratio = smallest/diagonal
			mpq_set_den(ratio,x[col]);
			mpq_abs(ratio,ratio);					        // ratio = |ratio|
			if ( mpq_cmp(ratio,tol) >= 0)				    // Is ratio >= tol?
				pivot = col;
			mpq_clear(tol); mpq_clear(ratio);			    // Free memory
		}
	}
	else if (order == 4)							        // Tolerance with largest pivot
	{
		pivot = SLIP_get_largest_pivot(x, pivs, n, top, xi);
		if ( mpz_sgn(x[col]) != 0 && pivs[col] < 0)			// Check x[col] vs largest potential pivot
		{
			mpq_t tol, ratio;
			mpq_inits(tol,ratio,NULL);
			mpq_set_d(tol, tolerance);				        // tol = user specified tolerance
			mpq_set_num(ratio,x[col]);				        // ratio = diagonal/largest
			mpq_set_den(ratio,x[pivot]);
			mpq_abs(ratio,ratio);					        // ratio = |ratio|
			if ( mpq_cmp(ratio,tol) >= 0)				    // Is ratio >= tol?
				pivot = col;
			mpq_clear(tol); mpq_clear(ratio);			    // Free memory
		}
	}
	else									                // Use the largest potential pivot
		pivot = SLIP_get_largest_pivot(x, pivs, n, top, xi);
	
	/* Reflect changes in row location & row_perm */
	if (pivot == -1) return pivot;						    // Error
	int intermed = pinv[pivot]; 						    // Must move pivot into position k
	int intermed2 = row_perm[k];
	
	/* Set row_perm[k] = pivot and row_perm[pinv[pivot]] = row_perm[k]
	   Also, set pinv[pivot] = k and pinv[row_perm[k]] = pinv[pivot] */
	row_perm[k] = pivot;
	row_perm[intermed] = intermed2;
	pinv[pivot] = k;
	pinv[intermed2] = intermed;
	pivs[pivot] = 1;							            // Row pivot is now pivotal
	mpz_set(rhos[k],x[pivot]); 						        // The kth pivot is x[pivot]
	return pivot;
}

/* Purpose: This function sorts the xi vector with respect to the current row permutation.
   This sort is efficient as its complexity is |x| log |x|.
   The idea of the sort is that you have xi[top, top+1, ...]. We essentially mask them
   and sort the masked vector (which is with respect to the row permutation). We then
   unmask them to get the correct value. For instance, the correct sorted order could
   be [4 2 5 1] because of the column permutation.
   Arguments:
   xi: nonzero pattern
   top: nonzeros are stored in xi[top..n-1]
   n: size of problem
   pinv: inverse row permutation
   row_perm: opposite of pinv. if pinv[j] = k then row_perm[k] = j */
void SLIP_sort_xi (int* xi, int top, int n, int* pinv, int* row_perm)
{
	for (int j = top; j < n; j++)		// Convert xi vector with respect to pinv
		xi[j] = pinv[xi[j]];
	std::sort(xi + top, xi + n);		// Sort xi[top..n-1]
	for (int j = top; j < n; j++)		// Place xi back in original value
		xi[j] = row_perm[xi[j]];
}

/* Purpose: This function performs the sparse REF triangular solve. i.e. (LD) x = A(:,k)
   The algorithm is described in the paper; however in essence it computes
   the nonzero pattern xi, then performs a sequence of IPGE operations on the
   nonzeros to obtain their final value. All operations are gauranteed to be integral.
   There are various enhancements in this code used to reduce the overall cost of the
   operations and minimize operations as much as possible.
   Arguments:
   L: partial L matrix
   A: input matrix
   k: iteration of algorithm
   xi: nonzero pattern vector
   q: column permutation
   rhos: sequence of pivots
   pinv: inverse row permutation
   row_perm: row permutation
   col_loc: column permutation
   h: history vector
   x: solution of system ==> kth column of L and U */
int SLIP_REF_triangular_solve(SLIP_mat* L, SLIP_mat* A, int k, int* xi, int* q,  mpz_t* rhos, int* pinv, int* row_perm, int* col_loc, int* h, mpz_t* x)
{
	int j, jnew, i, inew, p, m, top, n, col;
	/* Initialize REF TS by getting nonzero patern of x && obtaining A(:,k) */
	n = A->n;												                // Size of matrix and the dense vectors
	col = q[k];												                // Column we are solving for
	top = cs_reach(L, A, col, xi, pinv);									// Obtain nonzero pattern stored in xi[top..n]
	SLIP_sort_xi(xi, top, n, pinv, row_perm);								// Sort xi wrt sequence of pivots
	SLIP_reset_mpz_array_2(x,n,top,xi);									    // Reset x[i] = 0 for all i in nonzero pattern
	SLIP_reset_int_array_2(h,n,top,xi);									    // Reset h[i] = -1 for all i in nonzero pattern
	SLIP_get_column(A,col,x);										        // Set x = A(:,k)
	for (p = top; p < n; p++)										        // Iterate accross nonzeros in x
	{   
		/* Finalize x[j] */
		j = xi[p];											                // First nonzero term
		jnew = pinv[j];											            // Location of nonzero term
		if (mpz_sgn(x[j]) == 0) continue;								    // If x[j] == 0 no work must be done
		if (jnew < k)											            // jnew < k --> entries in U
		{
			/* History update */
			if (h[j] < jnew - 1)									        // HU must be performed
			{
				mpz_mul(x[j],x[j],rhos[jnew-1]);						    // x[j] = x[j] * rho[j-1]
				if (h[j] > -1)
					mpz_divexact(x[j],x[j],rhos[h[j]]);					    // x[j] = x[j] / rho[h[j]]
			}
			
			/* IPGE updates */
			col = col_loc[q[jnew]];									        // Corresponding column location of this x
			for (m = L->p[col]; m < L->p[col+1]; m++)						// Iterate accross nonzeros in Lij
			{
				i = L->i[m];									            // i value of Lij
				inew = pinv[i];									            // i location of Lij
				if (inew > jnew)
				{
					if (mpz_sgn(L->x[m]) == 0) continue;					// If lij is zero no update is performed
					if (mpz_sgn(x[i]) == 0)							        // lij is nonzero, x[i] is zero
					{
						/* x[i] = 0 --> Only perform IPGE update subtraction/division */
						if (jnew < 1)							            // No previous pivot
						{
							mpz_submul(x[i],L->x[m],x[j]);				    // x[i] = 0 - lij*x[j]
							h[i] = jnew;						            // Entry is up to date
						}
						else								                // Previous pivot exists
						{
							mpz_submul(x[i],L->x[m],x[j]);				    // x[i] = 0 - lij*x[j]
							mpz_divexact(x[i],x[i],rhos[jnew-1]);			// x[i] = x[i] / rho[j-1]
							h[i] = jnew;						            // Entry is up to date
						}
					}
					else									                // Both lij and x[i] are nonzero
					{
						/* x[i] != 0 --> History & IPGE update on x[i] */
						if (jnew < 1)							            // No previous pivot in this case
						{
							mpz_mul(x[i],x[i],rhos[0]);				        // x[i] = x[i]*rho[0]
							mpz_submul(x[i], L->x[m], x[j]);			    // x[i] = x[i] - lij*xj
							h[i] = jnew;						            // Entry is now up to date
						}
						else								                // There is a previous pivot
						{
							if (h[i] < jnew - 1)					        // History update if necessary
							{
								mpz_mul(x[i],x[i],rhos[jnew-1]);		    // x[i] = x[i] * rho[j-1]
								if (h[i] > -1)
									mpz_divexact(x[i],x[i],rhos[h[i]]);	    // x[i] = x[i] / rho[h[i]]
							}
							mpz_mul(x[i],x[i],rhos[jnew]);				    // x[i] = x[i] * rho[j]
							mpz_submul(x[i], L->x[m], x[j]);			    // x[i] = x[i] - lij*xj
							mpz_divexact(x[i],x[i],rhos[jnew-1]);			// x[i] = x[i] / rho[j-1] 
							h[i] = jnew;						            // Entry is up to date
						}
					}
				}
			}
		}
		else												                // Entries of L
		{
			/* History update */
			if (h[j] < k-1)
			{
				mpz_mul(x[j],x[j],rhos[k-1]);							    // x[j] = x[j] * rho[k-1]
				if (h[j] > -1)
					mpz_divexact(x[j],x[j],rhos[h[j]]);					    // x[j] = x[j] / rho[h[j]]
			}
		}
	}
	return top;												                // Output the beginning of nonzero pattern
}

/* Purpose: This function performs the SLIP LU factorization. This factorization is done via n iterations of the
   sparse REF triangular solve function. The overall factorization is PAQ = LDU
   Arguments:
   A: matrix to be factored
   L: lower triangular matrix
   U: upper triangular matrix
   S: stores guess on nnz and column permutation
   rhos: sequence of pivots
   pinv: inverse row permutation
   option: command options */
int SLIP_LU_Factor(SLIP_mat* A, SLIP_mat* L, SLIP_mat* U, SLIP_col* S,  mpz_t* rhos, int* pinv, SLIP_LU_Options* option)
{
	std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();	// Begin timing factorization
	int n, numRealloc = 0, k = 0, top, i, j, col, loc, lnz = 0, unz = 0, pivot, check, jnew, *xi, *h, *col_loc, *pivs, *row_perm;
	long size;

	/* Declare and initialize workspace */
	n = A->n;
	pivs = new int [n];					                    // Sequence of chosen pivots
	col_loc = new int [n];					                // Location of a column WRT the order
	h = new int [n];					                    // History vector
	xi = new int [2*n];					                    // Nonzero pattern
	row_perm = new int [n];					                // Row permutation, inverse of pinv
	SLIP_reset_int_array(pivs,n);
	SLIP_reset_int_array(h,n);
	
	/* Get most dense column and max of A */
	mpz_t sigma; mpz_init(sigma); mpz_set(sigma,A->x[0]);	// Initialize sigma
	for (i = 1; i < A->nz; i++)				                // Get sigma = max(A)
	{
		if(mpz_cmpabs(sigma,A->x[i]) < 0)
			mpz_set(sigma,A->x[i]);
	}
	mpz_abs(sigma,sigma);					                // sigma = |sigma|

	int gamma = A->p[1];
	for (i = 1; i<n; i++)					                // get gamma as most dense column
	{
		if( gamma < A->p[i+1] - A->p[i])
			gamma = A->p[i+1]-A->p[i];
	}
	mpfr_t temp; mpfr_init2(temp, 256); 
	mpfr_set_z(temp, sigma, MPFR_RNDN); 			        // temp = sigma
	
	/* Bound = gamma*log2(sigma sqrt(gamma)) */
	mpfr_mul_d(temp, temp, (double) sqrt(gamma), MPFR_RNDN);// temp = sigma*sqrt(gamma)
	mpfr_log2(temp, temp, MPFR_RNDN);			            // temp = log2(temp)
	double inner2 = mpfr_get_d(temp, MPFR_RNDN);		    // inner2 = temp
	mpfr_free_cache();					                    // Free cache from log2
	int bound = std::ceil(gamma*(inner2+1));		        // bound = gamma * inner2+1
	if (bound < 64) bound = 64;				                // Ensure bound is at least 64 bit
	mpfr_clear(temp); mpz_clear(sigma); 			        // Free memory
	
	/* Declare memory for x, L, and U */
	mpz_t* x = SLIP_initialize_mpz_array2(n,bound);	 		        // Initialize x
	for (i = 0; i < n; i++)					                // Initialize location based vectors
	{
		col_loc[S->q[i]] = i;
		pinv[i] = i;
		row_perm[i] = i;
	}
    SLIP_mat_alloc2(n, n, S->lnz, L);			            // Allocate L and U. 
    SLIP_mat_alloc2(n, n, S->unz, U);

	/* Iteration 0, must select pivot */
	col = S->q[0];
	check = 0;
	SLIP_get_column(A,col,x);				                // x = A(:,col)
	top = n-(A->p[col+1]-A->p[col]); j = 0;			        // top: nnz in column col
	for (i = A->p[col]; i < A->p[col+1]; i++)		        // Populate nonzero pattern
	{
		xi[top+j] = A->i[i];
		j+=1;
	}
	pivot = SLIP_get_pivot(x, pivs, n, top, xi, option->pivot, col, k, rhos, pinv, row_perm, option->tol);	// Get pivot
	if (pivot == -1) 					                    // The matrix is structurally deficient
	{	std::cout<<"\n**ERROR** No pivot at iteration 0! Matrix is singular!";
		return -1;
	}
	for (j = top; j < n; j++)				                // Populate L and U
	{
		jnew = xi[j];
		loc = pinv[jnew];				                    // Is j in U or L?
		if (loc <= k)					                    // U entries
		{
            U->i[unz] = jnew;                               // ith value of x[j]
			size = mpz_sizeinbase(x[jnew],2);	            // Allocate memory for x[j]
			mpz_init2(U->x[unz],size+2);                    // GMP manual: Allocated size should be size+2
			mpz_set(U->x[unz],x[jnew]);                     // Set U[x]
			unz += 1;                                       // Increment nnz of U
		}
		if (loc >= k)					                    // L entries
		{
			L->i[lnz] = jnew;                               // ith value of x[j]
			size = mpz_sizeinbase(x[jnew],2);
			mpz_init2(L->x[lnz],size+2);                    // GMP manual: Allocated size should be size+2
            mpz_set(L->x[lnz],x[jnew]);                     // Set L[x]
            lnz += 1;
		}
	}
	
	/* Iterations 1:n-1 (2:n in standard) */
	for (k = 1; k < n; k++)
	{
		L->p[k] = lnz;                                      // Column pointers for column k of L and U
        U->p[k] = unz;
		if (lnz + n > L->nzmax)                             // Reallocate memory if necessary
		{
            L->nz = lnz;                                    // Set L->nz = lnz
			SLIP_mat_realloc2(L);
			numRealloc +=1;
		}
		if (unz + n > U->nzmax)                             // Reallocate memory if necessary
		{
            U->nz = unz;                                    // Set U->nz = unz
			SLIP_mat_realloc2(U);
			numRealloc +=1;
		}
		top = SLIP_REF_triangular_solve(L, A, k, xi, S->q, rhos, pinv, row_perm, col_loc, h, x);                // LDx = A(:,k)
		pivot = SLIP_get_pivot(x, pivs, n, top, xi, option->pivot, S->q[k], k, rhos, pinv, row_perm, option->tol);	// Obtain pivot index
		if (pivot < 0) 							            // Error
		{
			std::cout<<"\n\n**ERROR** No pivot at iteration "<<k<<"! Matrix is singular!!\n\n";
			return -1;
		}
		for (j = top; j < n; j++) 					        // Iterate accross the nonzeros in x
		{
			jnew = xi[j];
			loc = pinv[jnew];					            // Location of x[j] in final matrix
			if (loc <= k)						            // loc <= k are rows above k, thus go to U
			{
				U->i[unz] = jnew;                           // Place the i location of the U->nz nonzero
				size = mpz_sizeinbase(x[jnew],2);
				mpz_init2(U->x[unz],size+2);                // GMP manual: Allocated size should be size+2
                mpz_set(U->x[unz],x[jnew]);                 // Place the x value of the U->nz nonzero
                unz += 1;                                   // Increment U->nz
			}
			if (loc >= k)						            // loc >= k are rows below k, thus go to L
			{
				L->i[lnz] = jnew;                           // Place the i location of the L->nz nonzero
				size = mpz_sizeinbase(x[jnew],2);
				mpz_init2(L->x[lnz],size+2);                // GMP manual: Allocated size should be size+2
                mpz_set(L->x[lnz],x[jnew]);                 // Place the x value of the L->nz nonzero
                lnz += 1;                                   // Increment L->nz
			}
		}
	}
	L->nz = lnz;
	U->nz = unz;
	L->p[n] = lnz;	    						            // Finalize L->p
	U->p[n] = unz;		    					            // Finalize U->p
	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	option->t_f = std::chrono::duration_cast<std::chrono::duration<float>>(t_end - t_begin);
	mpz_set(option->determinant, rhos[n-1]);				// Set the determinant of A (may be scaled)
	if (option->print == 1)							        // Print statistics if desired
	{
		std::cout<<"\n\n****Factorization Statistics****";
		std::cout<<"\n\nNum Nonzeros in A: "<<A->nz;
		std::cout<<"\nDimension of A, L, U: "<<n;
		std::cout<<"\nEstimated L nonzeros: "<<S->lnz;
		std::cout<<"\nActual L nonzeros: "<<L->nz;
		std::cout<<"\nEstimated U nonzeros: "<<S->unz;
		std::cout<<"\nActual U nonzeros: "<<U->nz;
		std::cout<<"\nNum Costly Reallocs: "<<numRealloc;
		std::cout<<"\n\n";
	}
	
	/* Free memory */
	SLIP_delete_mpz_array(x,n);
	delete[] xi; delete[] h; delete[] col_loc; delete[] pivs; delete[] row_perm;
	SLIP_mat_collapse2(L);							    // Collapse L 
	SLIP_mat_collapse2(U);							    // Collapse U
	return 0;								            // 0 indicates success
}

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------Sparse REF Forward and Backward Substitution------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* Purpose: This function performs sparse REF forward substitution with multiple right 
   hand side vectors. This is essentially the same as the sparse REF triangular 
   solve applied to each column of the right hand side vectors. Like the normal one,
   this function expects that the vector x is dense. As a result, the nonzero pattern 
   is not computed and each nonzero in x is iterated across. The system to solve is LDx = x
   Arguments:
   L: lower triangular matrix
   x: stores right hand side. Since it is multiple vectors it is stored as an mpz_t**
   rhos: sequence of pivots used in factorization
   numRHS: number of RHS vectors */
void SLIP_Forward_Sub(SLIP_mat* L, mpz_t** x, mpz_t* rhos, int numRHS)
{
	int i, j, p, k, n, m, mnew, **h;
	n = L->n;												                    // Size of x vector
	h = new int *[n];
	for (i = 0; i < n; i++)
	{
		h[i] = new int [numRHS];
		for (j = 0; j < numRHS; j++)
			h[i][j] = -1;
	}
	for (k = 0; k < numRHS; k++)										        // Iterate across each RHS vector
	{
		for (i = 0; i < n; i++)										            // Iterate accross all nonzeros in x. Assume x is dense
		{
			p = h[i][k];
			if (mpz_sgn(x[i][k]) == 0) continue;							    // Can skip these operations
			/* History Update */
			if (p < i-1)
			{
				mpz_mul(x[i][k],x[i][k],rhos[i-1]);						        // x[i] = x[i] * rhos[i-1]
				if (p > -1)						
					mpz_divexact(x[i][k],x[i][k],rhos[p]);					    // x[i] = x[i] / rhos[p]
			}
			/* IPGE updates */	
			for (m = L->p[i]; m < L->p[i+1]; m++)							    // Access the Lmi
			{
				mnew = L->i[m];									                // Location of Lmi
				if (mpz_sgn(L->x[m]) == 0) continue;						    // Can skip if Lx is zero
				if (mnew > i)									                // m > i
				{
					p = h[mnew][k];
					if (mpz_sgn(x[mnew][k]) == 0)						        // x[mnew] is zerp
					{
						mpz_submul(x[mnew][k], L->x[m], x[i][k]);				// x[m] = x[m] - lmi xi
						if (i > 0)
							mpz_divexact(x[mnew][k],x[mnew][k],rhos[i-1]);		// x[m] = x[m] / rhos[i-1]
					}
					else
					{
						if (p < i-1)							                // History update if necessary
						{
							mpz_mul(x[mnew][k],x[mnew][k],rhos[i-1]);		    // x[m] = x[m] * rhos[i-1]
							if (p > -1)
								mpz_divexact(x[mnew][k],x[mnew][k],rhos[p]);	// x[m] = x[m] / rhos[p]
						}
						mpz_mul(x[mnew][k],x[mnew][k],rhos[i]);				    // x[m] = x[m] * rhos[i]
						mpz_submul(x[mnew][k], L->x[m], x[i][k]);			    // x[m] = x[m] - lmi xi
						if (i > 0)
							mpz_divexact(x[mnew][k],x[mnew][k],rhos[i-1]);		// x[m] = x[m] / rhos[i-1]
					}
					h[mnew][k] = i;
				}
			}
		}
	}
	for (i = 0; i < n; i++)											            // Free h memory
		delete[] h[i];
	delete[] h;
}

/* Purpose: This function multiplies the x vector by the determinant of the matrix.
   Arguments:
   x: matrix to be multiplied
   det: constant to multiply by
   n: size of x
   numRHS: number of RHS vectors */
void SLIP_array_mul(mpz_t** x, mpz_t det, int n, int numRHS)
{
	for (int i = 0; i < n; i++)
		for (int k = 0; k < numRHS; k++)
			mpz_mul(x[i][k],x[i][k],det);
}

/* Purpose: This function divides the x vector by the determinant of the matrix. In general
   it can be used to divide an mpz array by an mpz constant. Note that the output
   is a rational mpq array. 
   Arguments:
   x2: solution of x/det
   x: input vector
   det: constant to divide by
   n: size of x and x2 
   numRHS: number of rhs vectors */
void SLIP_array_div(mpq_t** x2, mpz_t** x, mpz_t det, int n, int numRHS)
{
	mpq_t det2; mpq_init(det2); mpq_set_num(det2,det);	// Set det2 = det
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < numRHS; k++)
		{
			mpq_set_num(x2[i][k],x[i][k]);		        // Set x2[i] = x[i]
			mpq_div(x2[i][k],x2[i][k],det2);	        // x2[i] = x2[i] / det2
		}
	}
	mpq_clear(det2);					                // Free memory associated with det2
}

/* Purpose: This function performs sparse REF backward substitution. In essense it solves
   the sysem Ux = x. Note that prior to this, we expect x to be multiplied by
   the determinant of A. 
   Arguments:
   U: input upper triangular matrix
   x: right hand side vector
   numRHS: number of right hand side vectors */
void SLIP_Back_Sub (SLIP_mat* U, mpz_t** x, int numRHS)
{
	for (int k = 0; k < numRHS; k++)
	{
		for (int j = U->n-1; j >= 0; j--)				        // Start at x[n]
		{
			if (mpz_sgn(x[j][k]) == 0) continue;			    // If x[j] is zero stop
			mpz_divexact(x[j][k],x[j][k],U->x[U->p[j+1]-1]);	// Obtain x[j]
			for (int i = U->p[j]; i < U->p[j+1]-1; i++)
			{
				if (mpz_sgn(U->x[i])==0) continue;
				mpz_submul(x[U->i[i]][k], U->x[i], x[j][k]);	// x[i] = x[i] - U->x[i]*x[j]
			}
		}
	}
}

/* Purpose: This function solves the linear system LD^(-1)U x = b.
   Arguments:
   x: rational solution to the system
   b: right hand side vector
   rhos: sequence of pivots
   L: lower triangular matrix
   U: upper triangular matrix
   pinv: row permutation
   option: command options
   numRHS: number of RHS vectors */
void SLIP_Solve(mpq_t** x, mpz_t** b, mpz_t* rhos, SLIP_mat* L, SLIP_mat* U, int* pinv, SLIP_LU_Options* option, int numRHS)
{
	std::chrono::steady_clock::time_point t_s_begin = std::chrono::steady_clock::now();
	int i, k, n = L->n;
	mpz_t** b2 = SLIP_initialize_mpz_mat(n, numRHS);	    // Permuted b
	for (i = 0; i < L->nz; i++)				                // Permute entries in L
		L->i[i] = pinv[L->i[i]];
	for (i = 0; i < U->nz; i++)				                // Permute entries in U
		U->i[i] = pinv[U->i[i]];
	for (k = 0; k < numRHS; k++)				            // Set workspace b
	{
		for (i = 0; i < n; i++)			
			mpz_set(b2[pinv[i]][k],b[i][k]);
	}
	SLIP_Forward_Sub(L, b2, rhos, numRHS);		    // L*b2 = b2
	SLIP_array_mul(b2, rhos[n-1], n, numRHS);	    // b2 = b2 * det 
	SLIP_Back_Sub(U, b2, numRHS);			        // U b2 = b2
	SLIP_array_div(x, b2, rhos[n-1], n, numRHS);	// x = b2/det 
	std::chrono::steady_clock::time_point t_s_end = std::chrono::steady_clock::now();
	option->t_s = std::chrono::duration_cast<std::chrono::duration<float>>(t_s_end - t_s_begin);
	SLIP_delete_mpz_mat(b2, n, numRHS);			    // Free memory
}
#endif
