#ifndef SLIP_Output
#define SLIP_Output

/* This header includes functions for output */

/* This function prints the solution to the linear system to a file in either rational
   arithmetic, double precision, or variable precision floating point. The solution is
   permuted so that x is with respect to the original column permutation of A.
   Arguments:
   option: option struct which defines output data type
   x: the real solution vector
   n: the size of x */
void SLIP_LU_Print_File(SLIP_LU_Options* option, mpq_t** x, int n, int numRHS)
{
	std::string filename = option->out_file;		// Name of the output file
	std::ofstream output(filename);
	if (option->rat == 1)					        // Full precision rational arithmetic
	{
		output<<"Rational solution:\n";
		for (int i = 0; i < n; i++)
		{
			output <<"\n";
			for (int j = 0; j < numRHS; j++)
				output << x[i][j] << " ";
		}
	}
	else if (option->rat == 2)				        // Double precision 
	{
		double **x_doub = SLIP_initialize_double_mat(n, numRHS);
		for (int i = 0; i < n; i++)			        // Convert from rational to double
			for (int j = 0; j < numRHS; j++)
				x_doub[i][j] = mpq_get_d(x[i][j]);		
		output<<"Double solution:\n";
		for (int i = 0; i < n; i++)			        // Output the solution in double precision
		{
			output <<"\n";
			for (int j = 0; j < numRHS; j++)
				output << x_doub[i][j] << " ";
		}
		SLIP_delete_double_mat(x_doub, n, numRHS);
	}
	else							                // Variable precision floating point
	{
		mpfr_t** x_mpfr = SLIP_initialize_mpfr_mat(n, numRHS, option->prec);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < numRHS; j++)
				mpfr_set_q(x_mpfr[i][j], x[i][j], MPFR_RNDN);
		output << "Fixed precision solution of size " << option->prec << " bits\n";
		for (int i = 0; i < n; i++)
		{
			output << "\n";
			for (int j = 0; j < numRHS; j++)
			{
				char* outstring = NULL;				// Must do it this way...memory leaks otherwise
				mpfr_asprintf(&outstring, "%.*Rf", option->prec, x_mpfr[i][j]);		// Set outstring = x_mpfr[i] with prec digits of precision
				output << outstring << "\n";		// Output outstring
				mpfr_free_str(outstring);			// Free memory associated with outstring
			}
		}
		SLIP_delete_mpfr_mat(x_mpfr,n, numRHS);		// Free memory associated with x_mpfr
	}
	
}

/* This function prints statistics about the SLIP LU factorization
   Arguments:
   option: option struct telling how much info to print
   check: an int that tells whether the solution is correct or not
   nnz: number of nonzeros in L+U
   n: dimension of A
   numRHS: number of RHS vectors
   x: final solution vector */
void SLIP_LU_Print_Stats(SLIP_LU_Options* option, int check, int nnz, int n, int numRHS, mpq_t** x)
{
	if (option->check == 1)				// Was the solution checked?
	{
		if (check == 1)				    // 1 if correct
			std::cout<<"\nSolution is correct!";
		else					        // Incorrect solution. This should not happen
		{
			std::cout<<"\n****ERROR**** Solution incorrect!";
			std::cout<<"\n\nHave you modified the source code?";
			std::cout<<"\nReinstall and if this issue persists email me at: clouren@tamu.edu\n\n";
		}
		std::cout<<"\nSolution check time:\t\t"<<option->t_c.count();
	}
	option->t_tot = option->t_f + option->t_s+option->t_inp;	// Total Run time
	if (option->print2 == 1)			// Info about output file
	{
		std::cout<<"\nSolution output to file named:\t"<< option->out_file;
		if (option->rat == 1)
			std::cout<<"\nSolution output in full precision rational arithmetic";
		else if (option->rat == 2)
			std::cout<<"\nSolution output in double precision";
		else if (option->rat == 3)
			std::cout<<"\nSolution output in fixed precision of size: "<< option->prec << " bits";
		SLIP_LU_Print_File(option, x, n, numRHS);
	}
	if (option->print3 == 1)			// Print out factorization statistics
	{
        double nsquare = (double) n*n;
		std::cout<<"\nNumber of L+U nonzeros:\t\t"<<nnz;
		std::cout<<"\nL+U Density:\t\t\t" << (double) nnz/(nsquare);
		std::cout<<"\nInput manipulation time:\t" << option->t_inp.count();
		std::cout<<"\nSymbolic Analysis time:\t\t" << option->t_sym.count();
		std::cout<<"\nSLIP LU Factorization time:\t"<<option->t_f.count();
		std::cout<<"\nSLIP LU F/B Substitution time:\t"<<option->t_s.count();
		std::cout<<"\nSLIP LU Total Run Time:\t\t"<<option->t_tot.count()<<"\n\n";
	}
}
#endif