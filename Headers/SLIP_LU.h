#ifndef SLIP_LU
#define SLIP_LU

/* This is the primary header for SLIP LU. Include this in any code using SLIP LU
   Note: The other headers MUST be located in the same directory of SLIP_LU.h */
   
/* Purpose: This function will allow the user to take a matrix of their defined type (either double, mpfr_t, mpz_t, or mpq_t) and convert it 
   from their version of compressed column form to our data structure. The integrity of the user defined arrays are 
   maintained (therefore, one would need to delete these arrays)
   Arguments:
   p: The set of column pointers
   n: dimension of the matrix
   nz: number of nonzeros in the matrix (size of x and i vectors)
   i: set of row indices
   x_doub: Set of values in double precision. Pass NULL if not used
   x_int: Set of values stored as int. Pass NULL if not used
   x_mpz: Set of values in full precision integers. Pass NULL if not used
   x_mpq: Set of values as rational numbers. Pass NULL if not used
   prec: Precision associated with x_mpq. Pass any int if not used
   type: 0 mpz, 1 double, 2 int, 3 mpq, 4 mpfr  
   A: The input matrix. This structure should be allocated but not used yet 
   option: Used to compute the scaling parameter */	  
int SLIP_create_mat_ccf(int* p, int n, int nz, int* i, mpz_t* x_mpz, double *x_doub, int* x_int, mpq_t* x_mpq, 
		                mpfr_t* x_mpfr, int prec, int type, SLIP_mat* A, SLIP_LU_Options* option)
{
	if (type < 0 || type > 4)				// User input something wrong
	{
		std::cout<<"\nERROR: Type must be between 0-4 indicating the type of data structure used";
		return 0;
	}	
	else if (type == 0)					    // User has input a mpz array. All we must do is populate the structure
	{
		if (x_mpz == NULL)
		{
			std::cout<<"\nError! You selected mpz but passed nothing!";
			return 0;
		}
		SLIP_mpz_populate_mat(i, p, n, nz, x_mpz, A);
		mpq_set_ui(option->LU_scale, 1, 1);
		return 1;
	}
	else							        // User has input some type of array thats not easy
	{
		mpz_t* x_new = SLIP_initialize_mpz_array(nz);
		if (type == 1)					    // User has input a double array
		{
			if (x_doub == NULL)
			{
				std::cout<<"\nError! You selected double but passed nothing!";
				return 0;
			}
			SLIP_expand_double_array(x_doub, x_new, option->LU_scale, nz);
		}
		else if (type == 2)				    // User has input an int array
		{
			if (x_int == NULL)
			{
				std::cout<<"\nError! You selected int but passed nothing!";
				return 0;
			}
			SLIP_int_to_mpz(x_int, x_new, nz);
			mpq_set_ui(option->LU_scale,1,1);
		}
		else if (type == 3)			    	// User has input a mpq array
		{
			if (x_mpq == NULL)
			{
				std::cout<<"\nError! You selected mpq but passed nothing!";
				return 0;
			}
			SLIP_expand_mpq_array(x_mpq, x_new, option->LU_scale, nz);
		}
		else 						        // User has input an mpfr array
		{
			if (x_mpfr == NULL)
			{
				std::cout<<"\nError! You selected mpfr but passed nothing!";
				return 0;
			}
			SLIP_expand_mpfr_array(x_mpfr, x_new, option->LU_scale, nz, prec);
		}
		SLIP_mpz_populate_mat(i, p, n, nz, x_new, A);	// Create our matrix
		SLIP_delete_mpz_array(x_new, nz);		        // Free memory
		return 1;					                    // Success
	}
}

/* Purpose: This function will allow the user to take a matrix of their defined type (either double, mpfr_t, mpz_t, or mpq_t) and convert it from 
   their triplet form to our data structure. The integrity of the user defined arrays are 
   maintained (therefore, one would need to delete these arrays)
   Arguments:
   n: dimension of the matrix
   nz: number of nonzeros in the matrix (size of x, i, and j vectors)
   i: set of row indices
   j: set of column indices
   x_doub: Set of values in double precision. Pass NULL if not used
   x_int: Set of values stored as int. Pass NULL if not used
   x_mpz: Set of values in full precision integers. Pass NULL if not used
   x_mpq: Set of values as rational numbers. Pass NULL if not used
   prec: Precision associated with x_mpq. Pass any int if not used
   type: 0 mpz, 1 double, 2 int, 3 mpq, 4 mpfr  
   A: The input matrix. This structure should be allocated but not used yet
   option: Used for the scaling parameter (if necessary) */	  
int SLIP_create_mat_trip(int n, int nz, int* i, int* j, mpz_t* x_mpz, double *x_doub, int* x_int, mpq_t* x_mpq, 
		                 mpfr_t* x_mpfr, int prec, int type, SLIP_mat* A, SLIP_LU_Options* option)
{
	if (type < 0 || type > 4)				// User input something wrong
	{
		std::cout<<"\nERROR: Type must be between 0-4 indicating the type of data structure used";
		return 0;
	}	
	SLIP_trip *B;
	B = new SLIP_trip;
	B->n = n; B->m = n; B->nz = nz;	
	B->i = new int [nz];
	B->j = new int [nz];
	B->x = SLIP_initialize_mpz_array(nz);
	for (int k = 0; k < nz; k++)
	{
		B->i[k] = i[k];
		B->j[k] = j[k];
	}
	if (type == 0)					    // User has input a mpz array. All we must do is populate the structure
	{
		if (x_mpz == NULL)
		{
			std::cout<<"\nError! You selected mpz but passed nothing!";
			return 0;
		}
		for (int k = 0; k < nz; k++) mpz_set(B->x[k], x_mpz[k]);
		mpq_set_ui(option->LU_scale,1,1);
	}
	else if (type == 1)					// User has input a double array
	{	
		if (x_doub == NULL)
		{
			std::cout<<"\nError! You selected double but passed nothing!";
			return 0;
		}
		SLIP_expand_double_array(x_doub, B->x, option->LU_scale, nz);
	}
	else if (type == 2)				    // User has input an int array
	{
		if (x_int == NULL)
		{
			std::cout<<"\nError! You selected int but passed nothing!";
			return 0;
		}
		SLIP_int_to_mpz(x_int, B->x, nz);
		mpq_set_ui(option->LU_scale, 1,1);
	}
	else if (type == 3)				    // User has input a mpq array
	{
		if (x_mpq == NULL)
		{
			std::cout<<"\nError! You selected mpq but passed nothing!";
			return 0;
		}
		SLIP_expand_mpq_array(x_mpq, B->x, option->LU_scale, nz);
	}
	else 						        // User has input an mpfr array 
	{
		if (x_mpfr == NULL)
		{
			std::cout<<"\nError! You selected mpfr but passed nothing!";
			return 0;
		}
		SLIP_expand_mpfr_array(x_mpfr, B->x, option->LU_scale, nz, prec);
	}
	SLIP_trip_to_ccf(B,A); delete B;	// Create our matrix & clear memory
	return 1;
}

/* Purpose: This function will allow the user to take a rhs vector of their defined type (either double, mpfr_t, mpz_t, or mpq_t) and convert it from 
   their form to mpz. The integrity of the user defined arrays are maintained (therefore, one would need to delete these arrays)
   Arguments:
   n: dimension of the matrix/rhs
   b_doub: Set of values in double precision. Pass NULL if not used
   b_int: Set of values stored as int. Pass NULL if not used
   b_mpz: Set of values in full precision integers. Pass NULL if not used
   b_mpq: Set of values as rational numbers. Pass NULL if not used
   prec: Precision associated with x_mpq. Pass any int if not used
   type: 0 mpz, 1 double, 2 int, 3 mpq, 4 mpfr  
   b: The mpz_t data structure. This structure should be allocated but not used yet
   option: Used to scale the RHS if necessary */	  
int SLIP_create_rhs_2(int n, mpz_t* b_mpz, double *b_doub, int* b_int, mpq_t* b_mpq, mpfr_t* b_mpfr, int prec, int type, mpz_t* b, SLIP_LU_Options* option)
{
	if (type < 0 || type > 4)				// User input something wrong
	{
		std::cout<<"\nERROR: Type must be between 0-4 indicating the type of data structure used";
		return 0;
	}	
	if (type == 0)					        // User has input a mpz array. All we must do is populate the structure
	{
		if (b_mpz == NULL)
		{
			std::cout<<"\nError! You selected mpz but passed nothing!";
			return 0;
		}
		for (int k = 0; k < n; k++) mpz_set(b[k], b_mpz[k]);
		mpq_set_ui(option->b_scale,1,1);
	}
	else if (type == 1)				    	// User has input a double array
	{	
		if (b_doub == NULL)
		{
			std::cout<<"\nError! You selected double but passed nothing!";
			return 0;
		}
		SLIP_expand_double_array(b_doub, b, option->b_scale, n);
	}
	else if (type == 2)				        // User has input an int array
	{
		if (b_int == NULL)
		{
			std::cout<<"\nError! You selected int but passed nothing!";
			return 0;
		}
		SLIP_int_to_mpz(b_int, b, n);
		mpq_set_ui(option->b_scale, 1, 1);
	}
	else if (type == 3)				        // User has input a mpq array
	{
		if (b_mpq == NULL)
		{
			std::cout<<"\nError! You selected mpq but passed nothing!";
			return 0;
		}
		SLIP_expand_mpq_array(b_mpq, b, option->b_scale, n);
	}
	else 					            	// User has input an mpfr array
	{
		if (b_mpfr == NULL)
		{
			std::cout<<"\nError! You selected mpfr but passed nothing!";
			return 0;
		}
		SLIP_expand_mpfr_array(b_mpfr, b, option->b_scale, n, prec);
	}
	return 1;
}

/* Purpose: This function will allow the user to take a rhs matrix of their defined type (either double, mpfr_t, mpz_t, or mpq_t) and convert it from 
   their form to mpz. The integrity of the user defined arrays are maintained (therefore, one would need to delete these arrays). Note
   that this function assumes multiple RHS vectors
   Arguments:
   n: dimension of the matrix/rhs
   numRHS: number of RHS vectors
   b_doub: Set of values in double precision. Pass NULL if not used
   b_int: Set of values stored as int. Pass NULL if not used
   b_mpz: Set of values in full precision integers. Pass NULL if not used
   b_mpq: Set of values as rational numbers. Pass NULL if not used
   prec: Precision associated with x_mpq. Pass any int if not used
   type: 0 mpz, 1 double, 2 int, 3 mpq, 4 mpfr  
   b: The mpz_t data structure. This structure should be allocated but not used yet
   option: Used to scale the RHS if necessary */	  
int SLIP_create_rhs(int n, int numRHS, mpz_t** b_mpz, double **b_doub, int** b_int,  
			    mpq_t** b_mpq, mpfr_t** b_mpfr, int prec, int type, mpz_t** b, SLIP_LU_Options* option)
{
	if (type < 0 || type > 4)				// User input something wrong
	{
		std::cout<<"\nERROR: Type must be between 0-4 indicating the type of data structure used";
		return 0;
	}	
	if (type == 0)					        // User has input a mpz array. All we must do is populate the structure
	{
		if (b_mpz == NULL)
		{
			std::cout<<"\nError! You selected mpz but passed nothing!";
			return 0;
		}
		for (int i = 0; i < n; i++) 
			for (int j = 0; j < numRHS; j++)
				mpz_set(b[i][j], b_mpz[i][j]);
		mpq_set_ui(option->b_scale,1,1);
	}
	else if (type == 1)					    // User has input a double array
	{	
		if (b_doub == NULL)
		{
			std::cout<<"\nError! You selected double but passed nothing!";
			return 0;
		}
		SLIP_expand_double_mat(b_doub, b, option->b_scale, n, numRHS);
	}
	else if (type == 2)				        // User has input an int array
	{
		if (b_int == NULL)
		{
			std::cout<<"\nError! You selected int but passed nothing!";
			return 0;
		}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < numRHS; j++)
				mpz_set_si(b[i][j], b_int[i][j]);
		mpq_set_ui(option->b_scale, 1, 1);
	}
	else if (type == 3)				        // User has input a mpq array
	{
		if (b_mpq == NULL)
		{
			std::cout<<"\nError! You selected mpq but passed nothing!";
			return 0;
		}
		SLIP_expand_mpq_mat(b_mpq, b, option->b_scale, n, numRHS);
	}
	else 						            // User has input an mpfr array
	{
		if (b_mpfr == NULL)
		{
			std::cout<<"\nError! You selected mpfr but passed nothing!";
			return 0;
		}
		SLIP_expand_mpfr_mat(b_mpfr, b, option->b_scale, n, numRHS, prec);
	}
	return 1; 
}

/* Purpose: This function permutes the b vector & A matrix for the UMFPACK ordering.
   This is necessary since we use UMFPACK to preorder A, thus we must properly maintain b as well.
   Arguments:
   b: input right hand sides
   A: input matrix
   p_umf: row permutation given by UMFPACK 
   numRHS: number of RHS vectors */
void SLIP_UMFPACK_permute (mpz_t** b, SLIP_mat* A, int* p_umf, int numRHS)
{
	int k, i, index, n = A->n;
	int* pinv = new int [n];
	mpz_t** b2 = SLIP_initialize_mpz_mat(n, numRHS);
	for (k = 0; k < n; k++)
	{
		index = p_umf[k];
		pinv[index] = k;
	}
	for (k = 0; k < A->nz; k++)
		A->i[k] = pinv[A->i[k]];
	for (i = 0; i < numRHS; i++)
		for (k = 0; k < n; k++)
			mpz_set(b2[pinv[k]][i], b[k][i]);
	for (i = 0; i < numRHS; i++)
		for (k = 0; k < n; k++)
			mpz_set(b[k][i], b2[k][i]);
	SLIP_delete_mpz_mat(b2,n,numRHS);
	delete[] pinv;
}

/* Purpose: This function performs the symbolic ordering for SLIP LU. Currently, there 
   are four options: user defined order, COLAMD, AMD, or UMFPACK. 
   Arguments:
   A: Input matrix
   S: Structure which stores lnz, unz, and the column permutation
   option: Control parameters
   b: right hand side vectors (UMFPACK only)
   numRHS: Number of RHS vectors (UMFPACK only) */
void SLIP_LU_Symbolic(SLIP_mat* A, SLIP_col* S, SLIP_LU_Options* option, mpz_t** b, int numRHS) 
{
	std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
	int n = A->n, nz = A->nz, i, k;
	if (option->print == 1) SLIP_LU_Info();			    // Print info if needed
	
	/* No ordering */
	if (option->order == 2)
	{
		for (int i = 0; i < n+1; i++)
			S->q[i] = i;
		S->lnz = 10*nz; S->unz = S->lnz;	    	    // Guesses for number of L and U nonzeros
	}
	/* AMD */
	else if (option->order == 1)
	{
		double Control [AMD_CONTROL];		    	    // Declare AMD control
		amd_defaults(Control);				            // Set AMD defaults
		double Info [AMD_INFO];
		amd_order(n, A->p, A->i, S->q, Control, Info);	// Perform AMD
		S->lnz = Info[AMD_LNZ]; S->unz = S->lnz;	    // Guess for unz and lnz
		if (option->print == 1)				            // Output AMD info if desired
		{
			std::cout<<"\n****Column Ordering Information****\n";
			amd_control(Control);
			amd_info(Info);
		}
	}
	/* UMFPACK */
	else if (option->order == 3)
	{
		int status, sys, do_recip;
		option->pivot = 1;								// Must override pivoting scheme when using UMFPACK 
		int *p = new int [n];
		int *q = new int [n];
		double *Ax = new double [nz];					// Numeric values of Ax in double precision
		for (i = 0; i < nz; i++)
			Ax[i] = mpz_get_d(A->x[i]);
		double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
		void* Symbolic, *Numeric;
		/* Do UMFPACK */
		umfpack_di_defaults(Control);														// Set defaults
		Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;							// Unsymmetric ordering strategy
		Control[UMFPACK_PIVOT_TOLERANCE] = 0.000001;										// Pivot tolerance of 10^-6
		Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;									// Evaluates AMD, COLAMD, Metis and Nested Dissection
		status = umfpack_di_symbolic(n, n, A->p, A->i, Ax, &Symbolic, Control, Info);		// UMFPACK Symbolic
		status = umfpack_di_numeric(A->p, A->i, Ax, Symbolic, &Numeric, Control, Info);		// UMFPACK Numeric
		status = umfpack_di_get_numeric (NULL, NULL, NULL, NULL, NULL, NULL, p, q, NULL,
		&do_recip, NULL, Numeric);															// Get P and Q of UMFPACK
		umfpack_di_free_symbolic (&Symbolic) ;												// Free memory of symbolic
		umfpack_di_free_numeric (&Numeric) ;												// Free memory of numeric
	
		for (k = 0; k < n; k++)																// Set S->q as UMFPACK q
			S->q[k] = q[k];
		S->lnz = 2 * Info[UMFPACK_LNZ];														// UMFPACK lnz < SLIP lnz usually
		S->unz = 2 * Info[UMFPACK_UNZ];														// UMFPACK unz < SLIP unz usually
		SLIP_UMFPACK_permute(b, A, p, numRHS);
		delete[] Ax; delete[] p; delete[] q;
	}
	/* COLAMD */
	else
	{
		int Alen = 2*A->nz + 6 *(n+1) + 6*(n+1) + n;	// Declared as per COLAMD documentation
		int* A2 = new int [Alen];
		for (i = 0; i < n+1; i++)				        // Initialize S->q as per COLAMD documentation
			S->q[i] = A->p[i];
		for (i = 0; i < nz; i++)			 	       	// Initialize A2 per COLAMD documentation
			A2[i] = A->i[i];
		int stats [COLAMD_STATS];
		colamd(n, n, Alen, A2, S->q, (double *) NULL, stats);
		S->lnz = 10*A->nz; S->unz = S->lnz;		        // Guess for lnz and unz
		if (option->print == 1)			            	// Print stats if desired
		{
			std::cout<<"\n****Column Ordering Information****\n";
			colamd_report(stats);
			std::cout<<"\nEstimated L and U nonzeros: "<<S->lnz;
		}
		delete[] A2;
	}
	/* Make sure too much space isn't allocated */
	if (S->lnz > (double) n*n)					         // Guess exceeds max number of nnz in A
	{
		int nnz = ceil(0.5*n*n);
		S->lnz = nnz; S->unz = nnz;		
	}
	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	option->t_sym = std::chrono::duration_cast<std::chrono::duration<float>>(t_end - t_begin);
}

/* Purpose: This function returns the determinant of the input matrix A.
   Arguments:
   option: contains determinant of scaled A
   deter: output is mpq_t determinant
   n: size of matrix */
void SLIP_get_determinant(SLIP_LU_Options* option, mpq_t deter, int n)
{
	if (mpq_cmp_ui(option->LU_scale, 1, 1) == 0)	// Is LU scale 1?
		mpq_set_z(deter, option->determinant);
	else
	{
		mpq_set_z(deter, option->determinant);		// Set deter = deter of scaled A
		for (int k = 0; k < n; k++)			        // Scale deter
			mpq_div(deter, deter, option->LU_scale);
	}
}

/* Purpose: This function converts the SLIP LU L and U matrix into a Doolittle L and U
   matrix for output in compressed column form. 
   Arguments:
   pL: Doolittle L's p
   iL: Doolittle L's i
   xL: Doolittle L's x
   pU: Doolittle U's p
   iU: Doolittle U's i
   xU: Doolittle U's x
   L: SLIP's L
   U: SLIP's U
   rhos: Sequence of pivots
   option: command parameters */
void SLIP_get_Doolittle(int* pL, int* iL, mpq_t* xL, int* pU, int* iU, mpq_t* xU, SLIP_mat *L, SLIP_mat* U, mpz_t* rhos, SLIP_LU_Options* option)
{
	int n, lnz, unz, i, j, k;
	n = L->n; lnz = L->nz; unz = U->nz;
	mpq_t temp; mpq_init(temp); 
	
	/* As per Escobedo et. al. Column j of L is divided by rho[j] */
	for (k = 0; k < n; k++)
	{
		pL[k] = L->p[k];			               	// Populate pL and pU
		pU[k] = U->p[k];
		mpq_set_z(temp, rhos[k]);
		for (j = L->p[k]; j < L->p[k+1]; j++)
		{
			iL[j] = L->i[j];		              	// Populate iL
			mpq_set_z(xL[j], L->x[j]);
			mpq_div(xL[j], xL[j], temp);	     	// Populate xL
		}
	}
	pL[n] = L->p[n];				             	// p[n]
	pU[n] = U->p[n];
	/* As per Escobedo et. al., row j of U is divided by rho[j-1] */
	for (j = 0; j < unz; j++)
	{
		iU[j] = U->i[j];			            	// Populate iU
		mpq_set_z(xU[j], U->x[j]);
		if (iU[j] > 0)
		{
			mpq_set_z(temp, rhos[iU[j]-1]);			
			mpq_div(xU[j], xU[j], temp);
		}
		mpq_div(xU[j],xU[j],option->LU_scale);		// Populate xU
	}
	mpq_clear(temp);  
}


/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------Factor & Solve with Multiple RHS vectors----------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* Purpose: This function permutes x to get it back in its original form. That is x = Q*x.
   Arguments:
   x: Solution vector
   n: Size of solution vector
   numRHS: number of RHS vectors
   q: column permutation   */
void SLIP_Permute_x(mpq_t** x, int n, int numRHS, int* q)
{
	mpq_t** x2 = SLIP_initialize_mpq_mat(n, numRHS);	// Declare temp x
	for (int i = 0; i < n; i++)				            // Set x2 = Q*x
		for (int j = 0; j < numRHS; j++)
			mpq_set(x2[q[i]][j], x[i][j]);
	for (int i = 0; i < n; i++)				            // Return x = x2
		for (int j = 0; j < numRHS; j++)
			mpq_set(x[i][j],x2[i][j]);
	SLIP_delete_mpq_mat(x2,n,numRHS);
}

/* Purpose: This function scales the x matrix if necessary
   Arguments: 
   x: Solution matrix
   n: Size of A
   numRHS: number of RHS vectors
   option: contains scaling parameters */
void SLIP_Scale_x(mpq_t** x, int n, int numRHS, SLIP_LU_Options* option)
{
	if (option->LU_scale != NULL && (mpq_cmp_ui(option->LU_scale, 1, 1)) != 0)
		for (int i = 0; i < n; i++)
			for (int j = 0; j < numRHS; j++)
				mpq_mul(x[i][j], x[i][j], option->LU_scale);
	if (option->b_scale != NULL && (mpq_cmp_ui(option->b_scale, 1, 1)) != 0)
		for (int i = 0; i < n; i++)
			for (int j = 0; j < numRHS; j++)
				mpq_div(x[i][j], x[i][j], option->b_scale);
}

/* Purpose: This code utilizes the SLIP LU factorization. Soln is output as mpq_t mat
   Arguments:
   A: Compressed column form full precision matrix A
   S: Column ordering 
   b: Right hand side vectors stored as mpz_t matrix
   x: Solution vector stored as an mpq_t array
   numRHS: number of RHS vectors
   option: Control parameters */
int SLIP_LU_mpq(SLIP_mat* A, SLIP_col* S, mpz_t** b, mpq_t** x, int numRHS, SLIP_LU_Options* option)
{
	/* Declare memory */
	int check2, *pinv, n = A->n;
	SLIP_mat *L, *U;		
	L = new SLIP_mat; U = new SLIP_mat;		
	pinv = new int [n];
	mpz_t* rhos = SLIP_initialize_mpz_array(n);
	
	/* LU Factorization */
	check2 = SLIP_LU_Factor(A, L, U, S, rhos, pinv, option);
	if (check2 < 0)	return -1;
	/* FB Substitution */
	SLIP_Solve(x, b, rhos, L, U, pinv, option, numRHS);
	SLIP_Permute_x(x, n, numRHS, S->q);
	/* Check solution */
	if (option->check == 1)
		check2 = SLIP_LU_Check(A, x, n, numRHS, b, S->q, option);
	SLIP_Scale_x(x, n, numRHS, option);
	/* Print stats, free memory */
	SLIP_LU_Print_Stats(option, check2, L->nz+U->nz -n, n, numRHS, x);
	SLIP_delete_mpz_array(rhos,n);
	SLIP_mat_delete(L); SLIP_mat_delete(U);
	delete[] pinv;
	return 0;
}

/* This code utilizes the SLIP LU factorization. Soln is output as double matrix
   Arguments:
   A: Compressed column form full precision matrix A
   S: Column ordering 
   b: Right hand side vector stored as double array
   x: Solution vector stored as an mpq_t array
   numRHS: number of RHS vectors
   option: control parameters*/
int SLIP_LU_double(SLIP_mat* A, SLIP_col* S, mpz_t** b, double** x, int numRHS, SLIP_LU_Options* option)
{
	/* Declare memory */
	int check2, *pinv, n = A->n;
	mpq_t **x2 = SLIP_initialize_mpq_mat(n, numRHS);
	SLIP_mat *L, *U;	
	L = new SLIP_mat; U = new SLIP_mat;		
	pinv = new int [n];
	mpz_t* rhos = SLIP_initialize_mpz_array(n);
	
	/* LU factorization */
	check2 = SLIP_LU_Factor(A, L, U, S, rhos, pinv, option);
	if (check2 < 0)	return -1;
	/* FB Substitution */
	SLIP_Solve(x2, b, rhos, L, U, pinv, option, numRHS);
	SLIP_Permute_x(x2, n, numRHS, S->q);
	/* Check solution */
	if (option->check == 1)
		check2 = SLIP_LU_Check(A, x2, n, numRHS, b, S->q, option);
	SLIP_Scale_x(x2, n, numRHS, option);
	/* Output, free memory */
	for (int i = 0; i < n; i++)
		for (int j = 0; j < numRHS; j++)
			x[i][j] = mpq_get_d(x2[i][j]);
	SLIP_LU_Print_Stats(option, check2, L->nz+U->nz -n, n, numRHS, x2);
	SLIP_delete_mpq_mat(x2,n, numRHS);
	SLIP_delete_mpz_array(rhos,n);
	SLIP_mat_delete(L); SLIP_mat_delete(U);
	delete[] pinv;
	return 0; 
}

/* This code utilizes the SLIP LU factorization. Soln is output as mpfr_t
   Arguments:
   A: Compressed column form full precision matrix A
   S: Column ordering 
   b: Right hand side vector stored as mpz_t array
   x: Solution vector stored as an mpfr_t array
   numRHS: number of RHS vectors
   option: Control parameters */
int SLIP_LU_mpfr(SLIP_mat* A, SLIP_col* S, mpz_t** b, mpfr_t** x, int numRHS, SLIP_LU_Options* option)
{
	/* Declare memory */
	int check2, *pinv, n = A->n;
	mpq_t **x2 = SLIP_initialize_mpq_mat(n, numRHS);
	SLIP_mat *L, *U;
	L = new SLIP_mat; U = new SLIP_mat;		
	pinv = new int [n];
	mpz_t* rhos = SLIP_initialize_mpz_array(n);
	
	/* LU factorization */
	check2 = SLIP_LU_Factor(A, L, U, S, rhos, pinv, option);
	if (check2 < 0)	return -1;
	/* FB substituion */
	SLIP_Solve(x2, b, rhos, L, U, pinv, option, numRHS);
	SLIP_Permute_x(x2, n, numRHS, S->q);
	/* Check solution */
	if (option->check == 1)
		check2 = SLIP_LU_Check(A, x2, n, numRHS, b, S->q, option);
	SLIP_Scale_x(x2, n, numRHS, option);
	/* Output and free memory */
	for (int i = 0; i < n; i++)
		for (int j = 0; j < numRHS; j++)
			mpfr_set_q(x[i][j],x2[i][j], MPFR_RNDN);
	SLIP_LU_Print_Stats(option, check2, L->nz+U->nz -n, n, numRHS, x2);
	SLIP_delete_mpq_mat(x2,n, numRHS);
	SLIP_delete_mpz_array(rhos,n);
	SLIP_mat_delete(L); SLIP_mat_delete(U);
	delete[] pinv;
	return 0;
}

/* Purpose: This code utilizes the SLIP LU factorization. Soln is output as mpq_t vector
   This function is only intended for users who have a single RHS vector as input and output
   Arguments:
   A: Compressed column form full precision matrix A
   S: Column ordering 
   b: Right hand side vectors stored as mpz_t matrix
   x: Solution vector stored as an mpq_t array
   option: Control parameters */
int SLIP_LU_mpq2(SLIP_mat* A, SLIP_col* S, mpz_t* b, mpq_t* x, SLIP_LU_Options* option)
{
    int n = A->n;
    mpz_t ** b2 = SLIP_initialize_mpz_mat(n, 1);
    mpq_t ** x2 = SLIP_initialize_mpq_mat(n, 1);
    for (int k = 0; k < n; k++)
        mpz_set(b2[k][0], b[k]);
    int check = SLIP_LU_mpq(A, S, b2, x2, 1, option);
    for (int k = 0; k < n; k++)
        mpq_set(x[k], x2[k][0]);
    SLIP_delete_mpz_mat(b2,n,1);
    SLIP_delete_mpq_mat(x2,n,1);
    return check;
}

/* Purpose: This code utilizes the SLIP LU factorization. Soln is output as double vector
   This function is only intended for users who have a single RHS vector as input and output
   Arguments:
   A: Compressed column form full precision matrix A
   S: Column ordering 
   b: Right hand side vectors stored as mpz_t matrix
   x: Solution vector stored as an double array
   option: Control parameters */
int SLIP_LU_double2(SLIP_mat* A, SLIP_col* S, mpz_t* b, double* x, SLIP_LU_Options* option)
{
    int n = A->n;
    mpz_t ** b2 = SLIP_initialize_mpz_mat(n, 1);
    double ** x2 = SLIP_initialize_double_mat(n, 1);
    for (int k = 0; k < n; k++)
        mpz_set(b2[k][0], b[k]);
    int check = SLIP_LU_double(A, S, b2, x2, 1, option);
    for (int k = 0; k < n; k++)
        x[k] = x2[k][0];
    SLIP_delete_mpz_mat(b2,n,1);
    SLIP_delete_double_mat(x2,n,1);
    return check;
}

/* Purpose: This code utilizes the SLIP LU factorization. Soln is output as mpfr_t vector
   This function is only intended for users who have a single RHS vector as input and output
   Arguments:
   A: Compressed column form full precision matrix A
   S: Column ordering 
   b: Right hand side vectors stored as mpz_t matrix
   x: Solution vector stored as an double array
   option: Control parameters */
int SLIP_LU_mpfr2(SLIP_mat* A, SLIP_col* S, mpz_t* b, mpfr_t* x, SLIP_LU_Options* option)
{
    int n = A->n;
    mpz_t ** b2 = SLIP_initialize_mpz_mat(n, 1);
    mpfr_t ** x2 = SLIP_initialize_mpfr_mat(n, 1, option->prec);
    for (int k = 0; k < n; k++)
        mpz_set(b2[k][0], b[k]);
    int check = SLIP_LU_mpfr(A, S, b2, x2, 1, option);
    for (int k = 0; k < n; k++)
        mpfr_set(x[k], x2[k][0], MPFR_RNDN);
    SLIP_delete_mpz_mat(b2,n,1);
    SLIP_delete_mpfr_mat(x2,n,1);
    return check;
}

#endif