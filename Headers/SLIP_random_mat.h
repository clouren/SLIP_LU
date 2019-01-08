#ifndef SLIP_random
#define SLIP_random

/* This header contains routines for creating a dense, randomly generated matrix. The 
   purpose of this header is to compare SLIP LU with the dense REF LU factorization. 
   Note: The code to generate matrices here is exactly identical to the code used in REF LU*/

typedef struct random_mat_options
{
	int n;			    // Number of rows and columns
	int numRHS;		    // Number of RHS vectors
	int lower;		    // Integer lower bound
	int upper;		    // Integer upper bound
	int seed1;		    // Integer seed 1
	int seed2;		    // integer seed 2
	double density;		// Density of A
} random_mat_options;

/* Purpose: This function sets defaults to the random options struct
   Arguments:
   randopt: random options struct */
void SLIP_set_random_defaults(random_mat_options* randopt)
{
	randopt->n = DEFAULT_N;
	randopt->numRHS = DEFAULT_NUMRHS;
	randopt->lower = DEFAULT_LB;
	randopt->upper = DEFAULT_UB;
	randopt->seed1 = DEFAULT_SEED1;
	randopt->seed2 = DEFAULT_SEED2;
	randopt->density = DEFAULT_DENSITY;
}

/* Purpose: This function generates a random int matrix. This code is essentially the exact same as
   what is used in the software package "REF_LU" by Adolfo Escobedo which can be found at
   https://github.com/adolfoescobedo/REF-LU 
   Arguments:
   A: the dense int matrix to be populated
   randopt: options to generate A
   diag: a bool explaining how matrix is generated */   
void SLIP_generate_rand_mat(int** A, int n, int m, int lb, int ub, int seed1, int seed2, double density, bool diag)
{
	int i, j, index1, index2, temp, *row_index_rand;
	double frac_pos;							        // Fraction of entries expected to be positive
	int nnz = ceil(density*n*m);						// Number of nonzeros
	for (i = 0; i < m; i++)							    // Initialize A to all zeros
		for (j = 0; j < n; j++)
			A[i][j] = 0;
	
	/* Determine the ratio of positive to negative entries */
	if (lb >= 0) frac_pos = 1;
	else if (ub <= 0) frac_pos = 0;
	else frac_pos = (double) ub/(ub-lb);
	
	/* First, we ensure the matrix is nonsingular by placing a nonzero element along each column without repeating rows */
	srand(seed1);
	row_index_rand = new int [m];		
	for (i = 0; i < m; i++)							    // Populate row_index_rand with the natural order first
		row_index_rand[i] = i;
	
	if (!diag)								            // Generate the nonzeros to ensure nonsingularity RANDOMLY
		for (i = 0; i < 2*m; i++)
		{
			index1 = rand() % m;					    // Generate the random indices
			index2 = rand() % m;
			temp = row_index_rand[index1];				// Swap
			row_index_rand[index1] = row_index_rand[index2];
			row_index_rand[index2] = temp;
		}
	srand(2*seed1);								        // Generate new set of random numbers
	for (i = 0; i < std::min(n,m); i++)					// Populate the entries
	{
		index1 = row_index_rand[i];
		if ( (double) rand() / (RAND_MAX) < frac_pos)	// Check if this number should be positive or negative
			A[index1][i] = 1;
		else
			A[index1][i] = -1;
		nnz-=1;								            // Decrement nnz
	}
	delete[] row_index_rand;						    // Clear memory associated with row_index_rand

	/* Symbolically place the new values */
	while (nnz > 0)
	{
		index1 = rand() % m;						    // Randomly generate the indices
		index2 = rand() % n;
		if ( A[index1][index2] == 0)					// Place a symbolic entry
		{
			if ( (double) rand() / (RAND_MAX) < frac_pos)	// Check if number should be positive or negative
				A[index1][index2] = 1;
			else
				A[index1][index2] = -1;
			nnz -= 1;						            // Decrement nnz
		}
	}
	
	/* Generate the numerical entries */
	for (j = 0; j < n; j++)
		for (i = 0; i < m; i++)
		{ 
			if (A[i][j] < 0)					        // Populate with a negative value
				A[i][j] = -1 * ( rand() % (-1*lb)+1);
			else if (A[i][j] > 0)					    // Populate with positive value
				A[i][j] = rand() % ub + 1;
		}
}

/* Purpose: This function converts the dense int** A matrix into the SLIP_trip mat B
   Arguments:
   A: int mat to be covnerted
   B: SLIP_trip mat to be populated
   n: size of A */
void SLIP_dense_int_to_trip(int** A, SLIP_trip* B, int n)
{
	B->n = n; B->m = n;
	int nz = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (A[i][j] != 0)
				nz+=1;
	B->nz = nz;
	B->i = new int [nz];
	B->j = new int [nz];
	B->x = SLIP_initialize_mpz_array(nz);
	nz = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (A[i][j] != 0)
			{
				B->i[nz] = i;
				B->j[nz] = j;
				mpz_set_si(B->x[nz],A[i][j]);
				nz+=1;
			}
}

/* Purpose: This function processes the command line for SLIP random
   Arguments: 
   argc: Number of input arguments
   argv: Set of input parameters
   randopt: options */
int SLIP_process_rand_command(int argc, char* argv[], random_mat_options* randopt)
{
   for (int i = 1; i < argc; i++)
	{
		std::string arg = argv[i];
		if (arg == "n")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a number >=2 following n\n";
				return 0;
			}
			randopt->n = atoi(argv[i]);
			if (randopt->n < 2)
			{
				std::cout<<"\n****ERROR! Invalid n.\n\n";
				return 0;
			}
		}
		else if (arg == "b" || arg == "numRHS")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a number >=1 following b\n";
				return 0;
			}
			randopt->numRHS = atoi(argv[i]);
			if (randopt->numRHS < 1)
			{
				std::cout<<"\n****ERROR! Invalid b.\n\n";
				return 0;
			}
		}
		else if (arg == "l" || arg == "lower")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a number following l\n";
				return 0;
			}
			randopt->lower = atoi(argv[i]);
		}
		else if (arg == "u" || arg == "upper")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a number following u\n";
				return 0;
			}
			randopt->upper = atoi(argv[i]);
		}
		else if (arg == "s1" || arg == "seed1")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a number >=1 following s1\n";
				return 0;
			}
			randopt->seed1 = atoi(argv[i]);
			if (randopt->seed1 < 1)
			{
				std::cout<<"\n****ERROR! Invalid seed1.\n\n";
				return 0;
			}
		}
		else if (arg == "s2" || arg == "seed2")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a number >=1 following s2\n";
				return 0;
			}
			randopt->seed2 = atoi(argv[i]);
			if (randopt->seed2 < 1)
			{
				std::cout<<"\n****ERROR! Invalid seed2.\n\n";
				return 0;
			}
		}
		else if (arg == "d" || arg == "density")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a number (0,1] following d\n";
				return 0;
			}
			randopt->density = atof(argv[i]);
			if (randopt->density <= 0 || randopt->density > 1)
			{
				std::cout<<"\n****ERROR! Invalid density.\n\n";
				return 0;
			}
		}
		else
		{
			std::cout<<"\n\n**ERROR! Unknown command line parameter: "<<arg<<"\nIgnoring this parameter\n";
		}
	}
	return 1;
}

#endif