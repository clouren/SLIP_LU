#ifndef SLIP_Conversions
#define SLIP_Conversions

/* This header contains functions for conversions. Primarily converts matrices into
   compressed column form, arrays from mpq->mpfr, double->mpq, etc */

/* Purpose: p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c 
   This function comes from CSparse
   Arguments:
   p: vector to store the sum of c
   c: vector which is summed
   n: size of c */
int cs_cumsum (int *p, int *c, int n)
{
	int i, nz = 0 ;
	for (i = 0 ; i < n ; i++)
	{
		p [i] = nz ;
		nz += c [i] ;
		c [i] = p [i] ;
	}
	p [n] = nz ;
	return (nz) ;		    /* return sum (c [0..n-1]) */
}

/* Purpose: This function converts the triplet matrix B into the compressed column matrix A
   Arguments:
   B: matrix stored in triplet form
   A: matrix stored in ccf that will take B*/
void SLIP_trip_to_ccf(SLIP_trip* B, SLIP_mat* A)
{
	int k, p;
	SLIP_mat_alloc(B->n, B->m, B->nz, A);		// Allocate the sparse matrix A
	int* w = SLIP_initialize_int_array(A->n);	// Initialize w
	for (k = 0; k < B->nz; k++) w[B->j[k]]++;	// Column pointers
	cs_cumsum(A->p, w, A->n);			        // Column sums for A
	for (k = 0; k < B->nz; k++)
	{
		p = w[B->j[k]]++;
		A->i[p] = B->i[k];			            // Place values of i
		mpz_set(A->x[p],B->x[k]);		        // Place values of x
	}
	A->nz = A->nzmax;				            // Number of nonzeros in A
	SLIP_delete_mpz_array(B->x,B->nz);		    // Delete B and w
	delete[] B->i; delete[] B->j; delete[] w;
}

/* Purpose: This function converts a double array of size n to an appropriate mpz array
   of size n. To do this, the number is multiplied by 10^17 then, the GCD is found.
   This function allows the use of matrices in double precision to work with SLIP LU
   Arguments:
   x: double array that needs to be made integral
   x2: integral final array
   scale: the scaling factor used (x2 = scale * x)
   n: size of x */
void SLIP_expand_double_array(double* x, mpz_t* x2, mpq_t scale, int n)
{
	double expon = std::pow(10,17);				        // Double precision accurate ~17 decimals
	mpfr_t* x3 = SLIP_initialize_mpfr_array(n, 128);	// Quad precision in case input is huge
	mpq_set_d(scale,expon);					        	// scale = 10^17
	for (int i = 0; i < n; i++)
	{
		mpfr_set_d(x3[i], x[i], MPFR_RNDN);		    	// x3[i] = x[i]
		mpfr_mul_d(x3[i], x3[i], expon, MPFR_RNDN);		// x3[i] = x[i] * 10^17
		mpfr_get_z(x2[i], x3[i], MPFR_RNDN);			// x2[i] = x3[i]
	}
	SLIP_delete_mpfr_array(x3,n);				    	// Free memory associated with x3
		
	/* Compute the gcd to reduce the size of scale */
	mpz_t gcd, one;	mpz_inits(gcd, one, NULL);
	mpz_set(gcd, x2[0]); mpz_set_ui(one, 1);
	for (int i = 1; i < n; i++)				        	// Compute the GCD of the numbers
	{
		mpz_gcd(gcd,gcd,x2[i]);
		if ( (mpz_cmp(gcd,one) == 0))			    	// If gcd == 1 then stop
			break;
	}
	/* Scale all entries to make as small as possible */
	if ( (mpz_cmp(gcd,one) != 0))				    	// If gcd == 1 then stop
	{
		for (int i = 0; i < n; i++)
			mpz_divexact(x2[i],x2[i],gcd);
		mpq_t temp; mpq_init(temp);
		mpq_set_z(temp, gcd);
		mpq_div(scale, scale, temp);
		mpq_clear(temp);
	}
	mpz_clear(gcd); mpz_clear(one);				    	// Free memory 
}

/* Purpose: This function converts a mpfr array of size n and precision prec to an appropriate
   mpz array of size n. To do this, the number is multiplied by the appropriate power of 10 then the 
   gcd is found. This function allows mpfr arrays to be used within SLIP LU 
   Arguments:
   x: mpfr array to be expanded
   x2: full precision mpz array
   scale: scaling factor used (x2 = scale*x)
   n: size of x
   prec: associated precision used by mpfr */
void SLIP_expand_mpfr_array(mpfr_t* x, mpz_t* x2, mpq_t scale, int n, int prec)
{
	mpfr_t *x3 = SLIP_initialize_mpfr_array(n, prec);	// Create the new x array
	mpfr_t expon; mpfr_init2(expon, prec);
	mpfr_ui_pow_ui(expon, 10, prec, MPFR_RNDN);			// expon = 10^prec (overestimate)
	for (int i = 0; i < n; i++)
	{
		mpfr_mul(x3[i], x[i], expon, MPFR_RNDN);		// x3[i] = x[i]*10^prec
		mpfr_get_z(x2[i], x3[i], MPFR_RNDN);			// x2[i] = x3[i]
	}
	mpz_t temp_expon; mpz_init(temp_expon);
	mpfr_get_z(temp_expon, expon, MPFR_RNDN);
	mpq_set_z(scale, temp_expon);
	mpz_clear(temp_expon);
	SLIP_delete_mpfr_array(x3,n);				    	// Free memory
	mpfr_clear(expon);
	
	/* Find the gcd to reduce scale */
	mpz_t gcd, one; mpz_inits(gcd, one, NULL);
	mpz_set(gcd, x2[0]); mpz_set_ui(one, 1);
	for (int i = 1; i < n; i++)				        	// Compute the GCD of the numbers
	{
		mpz_gcd(gcd,gcd,x2[i]);
		if ( (mpz_cmp(gcd,one) == 0))			    	// If gcd == 1 stop
			break;
	}
	/* Scale all entries to make as small as possible */
	if ( (mpz_cmp(gcd,one) != 0))				    	// If gcd == 1 stop
	{
		for (int i = 0; i < n; i++)
			mpz_divexact(x2[i],x2[i],gcd);
		mpq_t temp; mpq_init(temp);
		mpq_set_z(temp,gcd);
		mpq_div(scale,scale,temp);
		mpq_clear(temp);
	}
	mpz_clear(gcd); mpz_clear(one);				    	// Free memory	
}

/* Purpose: This function converts a mpq array of size n into an appropriate mpz array of size n.
   To do this, the lcm of the denominators is found as a scaling factor. This function
   allows mpq arrays to be used in SLIP LU 
   Arguments:
   x: mpq array that needs to be converted
   x2: mpz array 
   scale: scaling factor. x2 = scale*x
   n: size of x */
void SLIP_expand_mpq_array(mpq_t* x, mpz_t* x2, mpq_t scale, int n)
{
	mpz_t temp; mpz_init(temp);
	mpz_t* x3 = SLIP_initialize_mpz_array(n);	// Initialize arrays
	mpq_t* x4 = SLIP_initialize_mpq_array(n);
	for (int i = 0; i < n; i++)			    	// x3 = denominators of x
		mpq_get_den(x3[i], x[i]);
	mpz_set(temp,x3[0]);
	for (int i = 1; i < n; i++)			    	// Find LCM of denominators of x
		mpz_lcm(temp, x3[i], temp);
	mpq_set_z(scale,temp);
	for (int i = 0; i < n; i++)	
	{
		mpq_mul(x4[i], x[i], scale);			// x4[i] = x[i]*temp
		mpz_set_q(x2[i], x4[i]);		    	// x2[i] = x4[i]
	}
	SLIP_delete_mpz_array(x3,n);				// Free memory
	SLIP_delete_mpq_array(x4,n);	
	mpz_clear(temp);
}

/* Purpose: This function populates the SLIP_mat A by the ccf stored vectors i, p, and x
   Arguments:
   i: row indices
   p: column pointers
   n: size of the matrix A
   nz: number of nonzeros in the matrix A
   x: set of values in A
   A: matrix to be populated */
void SLIP_mpz_populate_mat(int* i, int* p, int n, int nz, mpz_t* x, SLIP_mat* A)
{
	SLIP_mat_alloc(n, n, nz, A);		// Allocate the matrix A
	A->nz = nz;
	for (int k = 0; k <= n; k++)		// Set A->p
		A->p[k] = p[k];
	for (int k = 0; k < nz; k++)		// Set A->i and A->x
	{
		A->i[k] = i[k];
		mpz_set(A->x[k],x[k]);
	}
}

/* Purpose: This function converts the double array b1 into an mpz array
   Note that this will perform a truncation.
   Arguments:
   b1: double array
   b: mpz integral array
   n: size of b */
void SLIP_double_to_mpz_array (double* b1, mpz_t* b, int n)
{
	for (int i = 0; i < n; i++)
		mpz_set_d(b[i],b1[i]);
}

/* Purpose: This function converts mpq array to double
   This induces roundoff error via the final division
   Arguments:
   x2: double array
   x: mpq array
   n: size of b */
void SLIP_mpq_to_double (double* x2, mpq_t* x, int n)
{
	for (int i = 0; i < n; i++)
		x2[i] = mpq_get_d(x[i]);
}

/* Purpose: This function allows the user to obtain a solution using variable precision MPFR library
   Note, roundoff errors are induced in this final conversion
   Arguments:
   x: mpfr array
   x2: mpq array
   n: size of x */
void SLIP_mpq_to_mpfr (mpfr_t* x, mpq_t* x2, int n)
{
	for (int i = 0; i < n; i++)
		mpfr_set_q(x[i], x2[i], MPFR_RNDN);
}

/* Purpose: This function converts an int array to a mpz array
   Arguments:
   x: int array
   x2: mpz array
   n: dimension of x */
void SLIP_int_to_mpz (int* x, mpz_t* x2, int n)
{
	for (int i = 0; i < n; i++)
		mpz_set_si(x2[i],x[i]);
}

/* Purpose: This function converts a double matrix of size m*n to an appropriate mpz array
   of size m*n. To do this, the number is multiplied by 10^17 then, the GCD is found.
   This function allows the use of matrices in double precision to work with SLIP LU
   Arguments:
   x: double matrix that needs to be made integral
   x2: integral final matrix
   scale: the scaling factor used (x2 = scale * x)
   m: size of x
   n: size of x */
void SLIP_expand_double_mat(double** x, mpz_t** x2, mpq_t scale, int m, int n)
{
	int i, j;
	mpfr_t **x3 = SLIP_initialize_mpfr_mat(m, n, 128);	
	double expon = std::pow(10,17);								// expon = 10^17
	mpq_set_d(scale, expon);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			mpfr_set_d(x3[i][j], x[i][j], MPFR_RNDN);		    // x3[i][j] = x[i][j]
			mpfr_mul_d(x3[i][j], x3[i][j], expon, MPFR_RNDN);	// x3[i][j] = x[i][j]*10^17
			mpfr_get_z(x2[i][j], x3[i][j], MPFR_RNDN);		    // x2[i][j] = x3[i][j]
		}
	}
	SLIP_delete_mpfr_mat(x3, m, n);						        // Free memory of x3
	
	/* Compute the gcd to reduce the size of scale */
	mpz_t gcd, one;
	mpz_inits(gcd, one, NULL);
	mpz_set(gcd, x2[0][0]);
	mpz_set_ui(one, 1);
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)			                        // Compute the GCD of the numbers
		{
			mpz_gcd(gcd,gcd,x2[i][j]);
			if ( (mpz_cmp(gcd,one) == 0))	                    // If gcd == 1 then stop
				break;
		}
	/* Scale all entries to make as small as possible */
	if ( (mpz_cmp(gcd,one) != 0))			                    // If gcd == 1 stop
	{
		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++)
				mpz_divexact(x2[i][j],x2[i][j],gcd);
		mpq_t temp; mpq_init(temp);
		mpq_set_z(temp, gcd);
		mpq_div(scale, scale, temp);
		mpq_clear(temp);
	}
	mpz_clear(gcd); mpz_clear(one);			                    // Free memory
}

/* Purpose: This function converts a mpfr matrix of size m*n and precision prec to an appropriate
   mpz matrix of size m*n. To do this, the number is multiplied by the appropriate power of 10 then the 
   gcd is found. This function allows mpfr arrays to be used within SLIP LU 
   Arguments:
   x: mpfr matrix to be expanded
   x2: full precision mpz matrix
   scale: scaling factor used (x2 = scale*x)
   m: size of x
   n: size of x
   prec: associated precision used by mpfr */
void SLIP_expand_mpfr_mat(mpfr_t** x, mpz_t** x2, mpq_t scale, int m, int n, int prec)
{
	int i ,j;
	mpfr_t** x3 = SLIP_initialize_mpfr_mat(m, n, prec);
	mpfr_t expon; mpfr_init2(expon, prec);
	mpfr_ui_pow_ui(expon, 10, prec, MPFR_RNDN);			    // expon = 10^prec (overestimate)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			mpfr_mul(x3[i][j], x[i][j], expon, MPFR_RNDN);	// x3[i][j] = x[i][j]*expon
			mpfr_get_z(x2[i][j], x3[i][j], MPFR_RNDN);	    // x2[i][j] = x3[i][j]
		}
	}
	mpz_t temp_expon; mpz_init(temp_expon);
	mpfr_get_z(temp_expon, expon, MPFR_RNDN);
	mpq_set_z(scale, temp_expon);
	SLIP_delete_mpfr_mat(x3, m, n);
	mpfr_clear(expon);
	mpz_clear(temp_expon);
	
	/* Find the gcd to reduce scale */
	mpz_t gcd, one; mpz_inits(gcd, one, NULL);
	mpz_set(gcd, x2[0][0]);
	mpz_set_ui(one, 1);
	for (i = 0; i < m; i++)						            // Compute the GCD of the numbers
		for (j = 0; j < n; j++)
		{
			mpz_gcd(gcd,gcd,x2[i][j]);
			if ( (mpz_cmp(gcd,one) == 0))
				break;
		}
	/* Scale all entries to make as small as possible */
	if ( (mpz_cmp(gcd,one) != 0))					        // If gcd == 1 stop
	{
		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++)
				mpz_divexact(x2[i][j],x2[i][j],gcd);
		mpq_t temp; mpq_init(temp);
		mpq_set_z(temp, gcd);
		mpq_div(scale,scale,temp);
		mpq_clear(temp);
	}
	mpz_clear(gcd); mpz_clear(one);					        // Free memory	
}

/* Purpose: This function converts a mpq matrix of size m*n into an appropriate mpz array of size m*n.
   To do this, the lcm of the denominators is found as a scaling factor. This function
   allows mpq arrays to be used in SLIP LU 
   Arguments:
   x: mpq mat that needs to be converted
   x2: mpz mat
   scale: scaling factor. x2 = scale*x
   m: size of x
   n: size of x */
void SLIP_expand_mpq_mat(mpq_t** x, mpz_t** x2, mpq_t scale, int m, int n)
{
	int i , j;
	mpq_t** x4 = SLIP_initialize_mpq_mat(m, n);
	mpz_t** x3 = SLIP_initialize_mpz_mat(m, n);
	mpz_t temp; mpz_init(temp);
	for (i = 0; i < m; i++)					    // x3 = denominators of x
		for (j = 0; j < n; j++)
			mpq_get_den(x3[i][j], x[i][j]);
	mpz_set(temp,x3[0][0]);
	for (i = 0; i < m; i++)					    // Find LCM of denominators of x
		for (j = 0; j < n; j++)
			mpz_lcm(temp, x3[i][j], temp);
	mpq_set_z(scale,temp);
	for (i = 0; i < m; i++)					    // Make x integral
	{
		for (j = 0; j < n; j++)
		{
			mpq_mul(x4[i][j], x[i][j], scale);	// x4[i][j] = x[i][j]*temp
			mpz_set_q(x2[i][j], x4[i][j]);		// x2[i][j] = x4[i][j]
		}
	}
	
	SLIP_delete_mpz_mat(x3, m, n);
	SLIP_delete_mpq_mat(x4, m, n);
	mpz_clear(temp);
}

#endif
