#ifndef SLIP_mex 
#define SLIP_mex

/* This struct determines if the input matrix and right hand side vectors are 
   integral */
typedef struct int_input	// Is input matrix/RHS integral or not?
{
    double intA;			// > 0 if input A is integral
    double intb;			// > 0 if input b is integral
} int_input;

/* Purpose: This function checks the input of the matlab array. It ensures that
   the input matrix and right hand side vectors are correct dimension.
   Arguments:
   input: The matlab input array */
void SLIP_error_check2(const mxArray * input [])
{
    if (mxIsComplex (input[0]))		// Is the matrix real valued?
        mexErrMsgTxt ("Matrix must be real; try backslash instead") ;
    if (!mxIsSparse (input[0])) 	// Is the matrix sparse?
        mexErrMsgTxt ("Matrix must be sparse. Try again with A = sparse(A)") ;
    if (mxIsSparse (input[1])) 		// Is b sparse?
        mexErrMsgTxt ("Right hand side vector must be dense. Try again with b = full(b)") ;
    if (!mxIsStruct (input[2]))		// Is third argument the struct?
        mexErrMsgTxt ("Third argument must be the options struct") ;
    if (mxGetNumberOfFields(input[2]) != 9)
        mexErrMsgTxt("Error! The options struct must have 9 elements. Please reset it with option = SLIP_get_options;");
    int n = mxGetN(input[0]);
    if ( n != mxGetM(input[0]) )
        mexErrMsgTxt("A must be square");
    if ( n != mxGetM(input[1]) )
        mexErrMsgTxt("b is not of correct dimension");
}

/* Purpose: This function checks the input of the matlab array. It ensures that
   the input matrix is of correct dimension.
   Arguments:
   input: The matlab input array */
void SLIP_error_check3(const mxArray * input [])
{
	if (mxIsComplex (input[0]))		// Is the matrix real valued?
		mexErrMsgTxt ("Matrix must be real; try backslash instead") ;
    if (!mxIsSparse (input[0])) 		// Is the matrix sparse?
		mexErrMsgTxt ("Matrix must be sparse. Try again with A = sparse(A)") ;
	if (!mxIsStruct (input[1]))		// Is second argument the struct?
		mexErrMsgTxt ("Third argument must be the options struct") ;
	if (mxGetNumberOfFields(input[1]) != 9)
		mexErrMsgTxt("Error! The options struct must have 9 elements. Please reset it with option = SLIP_get_options;");
	if ( mxGetN (input[0]) != mxGetM (input[0]))				// Is the input square?
		mexErrMsgTxt("A must be square");
}


/* Purpose: This function reads in the necessary information from the options struct
   Arguments:
   input: The input A matrix, b vector and options 
   inter: determines if the input matrix is integer or not
   option: Control parameters */
void SLIP_get_matlab_options(const mxArray* input[], int_input* inter, SLIP_LU_Options* option)
{
	mxArray* tmp;
	
	tmp = mxGetFieldByNumber(input[2], 0, 0);		// Get the column ordering
	if (tmp == NULL)
		mexErrMsgTxt("Error for column ordering");
	double order = mxGetScalar(tmp);

	tmp = mxGetFieldByNumber(input[2], 0, 1);		// Get the pivoting scheme
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting pivot");
	double piv = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[2], 0, 2);		// Determine how much info is printed
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting printing parameter");
	option->print = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[2], 0, 3);		// Determine how much info is printed
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting printing parameter");
	option->print2 = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[2], 0, 4);		// Determine how much info is printed
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting printing parameter");
	option->print3 = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[2], 0, 5);		// Determine if final solution vector will be checked for accuracy
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting the check parameter");
	double checker = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[2], 0, 6);		// Set equal to 1 if input is integral
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting the int parameter");
	inter->intA = mxGetScalar(tmp);

	tmp = mxGetFieldByNumber(input[2], 0, 7);		// Set equal to 1 if input is integral
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting the intb parameter");
	inter->intb = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[2], 0, 8);		// Tolerance if some form of tolerance partial pivoting is used
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting the tolerance parameter");
	option->tol = mxGetScalar(tmp);
	
	/* Verify that the parameters are correct */
	if (order <= 3 && order >= 0) option->order = (int) order;
	if (piv <= 5 && piv >= 0) option->pivot = (int) piv;
	if ( checker > 0) option->check = 1;
	if (option->tol > 1 || option->tol <= 0) option->tol = DEFAULT_TOL;
}

/* Purpose: This function reads in the necessary information from the options struct
   Arguments:
   input: The input A matrix and options 
   option: Control parameters */
void SLIP_get_matlab_options3(const mxArray* input[], SLIP_LU_Options* option)
{
	mxArray* tmp;
	
	tmp = mxGetFieldByNumber(input[1], 0, 0);		// Get the column ordering
	if (tmp == NULL)
		mexErrMsgTxt("Error for column ordering");
	double order = mxGetScalar(tmp);

	tmp = mxGetFieldByNumber(input[1], 0, 1);		// Get the pivoting scheme
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting pivot");
	double piv = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[1], 0, 2);		// Determine how much info is printed
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting printing parameter");
	option->print = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[1], 0, 3);		// Determine how much info is printed
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting printing parameter");
	option->print2 = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[1], 0, 4);		// Determine how much info is printed
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting printing parameter");
	option->print3 = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[1], 0, 5);		// Determine if final solution vector will be checked for accuracy
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting the check parameter");
	double checker = mxGetScalar(tmp);
	
	tmp = mxGetFieldByNumber(input[1], 0, 8);		// Tolerance if some form of tolerance partial pivoting is used
	if (tmp == NULL)
		mexErrMsgTxt("Error at getting the tolerance parameter");
	option->tol = mxGetScalar(tmp);
	
	/* Verify that the parameters are correct */
	if (order <= 2 && order >= 0) option->order = (int) order;
	if (piv <= 5 && piv >= 0) option->pivot = (int) piv;
	if ( checker > 0) option->check = 1;
	if (option->tol > 1 || option->tol <= 0) option->tol = DEFAULT_TOL;
}

/* Purpose: This function converts an int to an INT 
   Arguments:
   x: the int vector
   y: the INT vector
   n: the size of x and y */
void SLIP_int_to_INT(int* x, INT* y, int n)
{
	for (int i = 0; i < n; i++)
		y[i] = (INT) x[i];
}

/* Purpose: This function converts an INT to an int
   Arguments:
   x: the int vector
   y: the INT vector
   n: the size of x and y */
void SLIP_INT_to_int(int* x, INT* y, int n)
{
	for (int i = 0; i < n; i++)
		x[i] = (int) y[i];
}

/* Purpose: This function checks if the input matrix & RHS has numbers too large for double
   Arguments:
   x: The numeric values of A 
   b: the numeric values of rhs vector
   nz: number of nonzeros in A
   n: dimension of A */
void SLIP_mex_check_for_inf(double* x, double* b, int nz, int n)
{
    for (int k = 0; k < nz; k++)
        if (mxIsInf(x[k]))
            mexErrMsgTxt("Numbers are too large for double. Please try the C++ code with mpfr, mpq, or mpz");
			
	for (int k = 0; k < n; k++)
		if (mxIsInf(b[k]))
			mexErrMsgTxt("Numbers are too large for double. Please try the C++ code with mpfr, mpq, or mpz");
}

/* Purpose: This function reads in the A matrix and right hand side vectors. 
   Arguments:
   A: Internal SLIP Mat stored in ccf 
   Ax: Numeric values of A from matlab
   Ap: Column pointers of A from matlab 
   Ai: row indices of A from matlab
   bx: multiple RHS vectors from matlab
   b2: mpz matrix used internally (b2 = bx)
   nA: number of columns in A
   nz: number of nonzeros in A
   nb: number of columns in b
   inter: if b and/or are already integral
   option: command options
   flag: determines if b exists or not */
void SLIP_mex_get_A_and_b(SLIP_mat* A, double* Ax, int* Ap, int* Ai, double** bx, mpz_t** b2, int nA, int nz, int nb, int_input* inter, SLIP_LU_Options* option, int flag)
{
    std::chrono::steady_clock::time_point t_inp_begin = std::chrono::steady_clock::now();   // Start time of input manipulation
    int test;
    /* Read in A */
    if (inter->intA > 0)                                    // Input is already integral
    {
        int* real_x = new int [nz];                              // Declare memory
        for (int k = 0; k < nz; k++) real_x[k] = (int) Ax[k];               
        test = SLIP_create_mat_ccf(Ap, nA, nz, Ai, NULL, NULL, real_x, NULL, NULL, 1, 2, A, option);         // Create A with no scaling
        delete[] real_x;
    }
    else
        test = SLIP_create_mat_ccf(Ap, nA, nz, Ai, NULL, Ax, NULL, NULL, NULL, 1, 1, A, option);             // Create A with scaling
    if (test == 0) 
        mexErrMsgTxt("Issue reading in A. Please ensure matrix is correct and try again");
    /* Read in b */
    if (flag != -1)
    {
        if (inter->intb > 0)                                                                                // Is b integral?
        {
            int** real_b = SLIP_initialize_int_mat(nA, nb);
            for (int k = 0; k < nb; k++)
                for (int j = 0; j < nA; j++)
                    real_b[j][k] = (int) bx[j][k]; 
            test = SLIP_create_rhs(nA, nb, NULL, NULL, real_b, NULL, NULL, 1, 2, b2, option);        // Create b
            SLIP_delete_int_mat(real_b, nA, nb);
        }
        else
            test = SLIP_create_rhs(nA, nb, NULL, bx, NULL, NULL, NULL, 1, 1, b2, option);            // Create b
        if (test == 0) 
            mexErrMsgTxt("Issue reading in b. Please ensure RHS is correct and try again");
    }
    delete[] Ap; delete[] Ai;
    std::chrono::steady_clock::time_point t_inp_end = std::chrono::steady_clock::now();
    option->t_inp = std::chrono::duration_cast<std::chrono::duration<float>>(t_inp_end - t_inp_begin);
}

/* Purpose: This function outputs a double array as a mxArray
   Arguments:
   x: The array to be output
   n: the size of x */
mxArray* SLIP_mex_output_double_array(double* x, int n)
{
	mxArray* Xmatlab = mxCreateDoubleMatrix (n, 1, mxREAL); 	// Create a n * 1 array
	double* x2 = mxGetPr(Xmatlab);					            // Get the numeric values
	for (int k = 0; k < n; k++) x2[k] = x[k];			        // Set Xmatlab = x
	return Xmatlab;
}

/* Purpose: This function outputs a double matrix as a mxArray
   Arguments:
   x: The matrix to be output
   m: size of x
   n: the size of x */
mxArray* SLIP_mex_output_double_mat(double** x, int m, int n)
{
	mxArray* Xmatlab = mxCreateDoubleMatrix (m, n, mxREAL); 	// Create a m * n array
	double* x2 = mxGetPr(Xmatlab);					            // Get the numeric values
	int count = 0;
	for (int j = 0; j< n; j++)                                  // Populate the nonzeros in output matrix
		for (int i = 0; i < m; i++)
		{
			x2[count] = x[i][j];
			count+=1;
		}
	return Xmatlab;
}

/* Purpose: This function outputs an int array as a mxArray. 
   Arguments:
   x: int array to be output
   n: size of x */
mxArray* SLIP_mex_output_int_array(int* x, int n)
{
	mxArray* Xmatlab = mxCreateDoubleMatrix (n, 1, mxREAL);	// Create a n*1 matlab array
	double* x2 = mxGetPr(Xmatlab);				        	// Numeric values of Xmatlab
	for (int k = 0; k < n; k++)
		x2[k] = (double) x[k];	    		            	// Get cast x as double
	return Xmatlab;
}

/* Purpose: This function outputs the p matrix from pinv as a mxArray 
   Arguments:
   pinv: pinv
   n: size of pinv */
mxArray* SLIP_mex_output_p(int* pinv, int n)
{
	mxArray* Pmatlab = mxCreateDoubleMatrix (n, 1, mxREAL);		// Create a n*1 array
	double* x = mxGetPr(Pmatlab);					// Numeric values of Pmatlab
	for (int k = 0; k < n; k++)
		x[pinv[k]] = (double) k;					// Set Pmatlab[pinv[k]] = k
	return Pmatlab;
}

/* Purpose: This function outputs a sparse matrix
   Arguments:
   n: Columns of matrix
   m: rows of matrix
   nz: number of nonzeros in matrix
   p: column pointers
   i: row indices
   x: numeric values */
mxArray* SLIP_mex_output_matrix(int n, int m, int nz, int* p, int* i, double* x)
{
	mxArray *Amatlab = mxCreateSparse( (INT) m, (INT) n, (INT) nz, mxREAL);	// Create a m*n sparse matrix
	INT* pA = new INT [n+1];
	INT* iA = new INT [nz];
	double* xA = new double [nz];
	pA = (INT*) mxGetJc(Amatlab);		// pA = A->p
	iA = (INT*) mxGetIr(Amatlab);		// iA = A->i
	xA = (double*) mxGetPr(Amatlab);	// xA = A->x
	
	for (int k = 0; k < n+1; k++)		// Populate A->p
		pA[k] = (INT) p[k];
	for (int k = 0; k < nz; k++)		// Populate A->x and A->i
	{
		iA[k] = (INT) i[k];
		xA[k] = x[k];
	}
	return Amatlab;
}
#endif