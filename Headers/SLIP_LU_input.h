#ifndef SLIP_input
#define SLIP_input 

/* This header is for command line arguments, reading in matrices, etc */

/*----------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
-----------------------------Command Line Arguments---------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------*/

/* Purpose: Set default parameters. These parameters are defined in SLIP_LU_config.h
   Arguments:
   option: The struct that will have defaults set */
void SLIP_Set_Options_Defaults(SLIP_LU_Options* option)
{
	option->check = DEFAULT_CHECK;
	option->pivot = DEFAULT_PIVOT;
	option->order = DEFAULT_ORDER;
	option->print = DEFAULT_PRINT;
	option->print2 = DEFAULT_PRINT2;
	option->print3 = DEFAULT_PRINT3;
	option->tol = DEFAULT_TOL;
	option->mat_name = DEFAULT_MAT;
	option->rhs_name = DEFAULT_RHS;
	option->out_file = DEFAULT_OUTFILE;
	option->rat = DEFAULT_RAT;
	option->prec = DEFAULT_PRECISION;
	std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();
	option->t_inp = std::chrono::duration_cast<std::chrono::duration<float>>(t - t);
	mpq_init(option->LU_scale); 		// Initialize the mpz elements
	mpq_init(option->b_scale);
	mpz_init(option->determinant);
	mpq_set_ui(option->LU_scale, 1,1);	// Initial scales are 1
	mpq_set_ui(option->b_scale, 1,1);
}

/* Purpose: This function shows the usage of the code. 
   Arguments: None*/
void SLIP_show_usage()
{
	std::cout<<"\n\n\t\t****USAGE****";
	std::cout<<"\n\t./SLIP_LU.exe followed by:";
	std::cout<<"\n\tc: Indicates soln will be checked\n\tp or piv: indicate type of pivoting";
	std::cout<<"\n\tcol or q: column order used: 0: colamd, 1: amd, 2: none";
	std::cout<<"\n\tt or tol: tolerance parameter \n\to2 or out2: output printed to screen";
	std::cout<<"\n\tf or file: filenames. must be of format MATRIX_NAME RHS_NAME";
	std::cout<<"\n\tof or outfile: output filename";
	std::cout<<"\n\to or out: output will be printed to file. Must be followed by 1: rational, 2: double, 3 PREC: float of precision PREC" ;
	std::cout<<"\n****REFER TO README.txt FOR DETAILED DESCRIPTION OF INPUT PARAMETERS****\n";
}

/* Purpose: This function prints out the user specified/default options 
   Arguments:
   option: struct containing all of the options */
void SLIP_print_options(SLIP_LU_Options* option)
{
	std::string piv, order;
	if (option->order == 0)
		order = "the COLAMD";
	else if (option->order == 1)
		order = "the AMD";
	else 
		order = "No";
	if (option->pivot == 0)
		piv = "smallest";
	else if (option->pivot == 1)
		piv = "diagonal";
	else if (option->pivot == 2)
		piv = "first nonzero";
	else if (option->pivot == 3)
		piv = "diagonal small tolerance";
	else if (option->pivot == 4)
		piv = "diagonal large tolerance";
	else
		piv = "largest";
	
	std::cout<<"\n\n****COMMAND PARAMETERS****";
	std::cout<<"\n\nSolving the linear system represented by the matrix "<<option->mat_name <<" and RHS "<<option->rhs_name;
	std::cout<<"\nUsing "<<order<<" ordering and selecting the "<<piv<<" pivot";
	if (option->pivot == 3 || option->pivot == 4)
		std::cout<<"\nTolerance used: "<<option->tol<<"\n";
	std::cout<<"\nOutput file is: "<<option->out_file;
}

/* Purpose: This processes the command line for user specified options
   Arguments:
   argc: number of command line arguments
   argv: set of command line arguments
   option: struct containing the command options */
int SLIP_LU_process_command_line(int argc, char* argv[], SLIP_LU_Options* option)
{
	for (int i = 1; i < argc; i++)
	{
		std::string arg = argv[i];
		if (arg == "help")
		{
			SLIP_show_usage();
			return 0;
		}
		else if (arg == "c" || arg == "check")
			option->check = 1;
		else if (arg == "p" || arg == "piv")
		{
			if (!argv[++i]) 
			{
				std::cout<<"\n****ERROR! There must be a pivot argument between 0-5 following p\n";
				return 0;
			}
			option->pivot = atoi(argv[i]);
			if (option->pivot < 0 || option->pivot > 5)
			{
				std::cout<<"\n****ERROR! Invalid pivot selection!\nDefaulting to smallest pivot\n\n";
				option->pivot = 0;
			}
		}
		else if (arg == "q" || arg == "col")    
		{
			if (!argv[++i])
			{
				std::cout<<"\n****ERROR! There must be an argument between 0-2 following q\n";
				return 0;
			}
			option->order = atoi(argv[i]);
			if (option->order < 0 || option->order > 2)
			{
				std::cout<<"\n****ERROR! Invalid column ordering\nDefaulting to COLAMD\n\n";
				option->order = 0;
			}
		}
		else if (arg == "t" || arg == "tol")
		{
			if (!argv[++i])
			{
				std::cout<<"\n****ERROR! There must be a non negative tolerance value following t\n";
				return 0;
			}
			else if (!atoi(argv[i]))
			{
				std::cout<<"\n****ERROR! There must be a non negative tolerance value following t\n";
				return 0;
			}
			option->tol = atof(argv[i]);
			if (option->tol < 0)
			{
				std::cout<<"\n****ERROR! Invalid Tolerance, tolerance must be non-negative\n";
				return 0;
			}
		}
		else if (arg == "out2" || arg == "o2")
			option->print = 1;
		else if (arg == "out" || arg == "o")
		{
			option->print2 = 1;
			if (!argv[++i])
			{
				std::cout<<"\n****ERROR! o or out must be followed by 1 (rational) 2 (double) or 3 (variable precision) \n";
				return 0;
			}
			else if (!atoi(argv[i]))
			{
				std::cout<<"\n****ERROR! o or out must be followed by 1 (rational) 2 (double) or 3 (variable precision) \n";
				return 0;
			}
			option->rat = atoi(argv[i]);
			if (option->rat < 1 || option->rat > 3)
			{
				std::cout<<"\n\n****ERROR! Invalid output type!\n Defaulting to rational";
				option->rat = 1;
			}
			if (option->rat == 3)
			{
				if (!argv[++i])
				{
					std::cout<<"\n****ERROR! Precision must be specified\n";
					return 0;
				}
				else if (!atoi(argv[i]))
				{
					std::cout<<"\n****ERROR! Precision must be specified\n";
					return 0;
				}
				option->prec = atoi(argv[i]);
				if (option->prec < 2)
				{
					std::cout<<"\n\n****ERROR! Invalid precision. p >= 2\n";
					return 0;
				}
			}
		}
		else if (arg == "f" || arg == "file")
		{
			if (!argv[++i])
			{
				std::cout<<"\n****ERROR! Matrix name must be entered\n";
				return 0;
			}
			option->mat_name = argv[i];
			if (!argv[++i])
			{
				std::cout<<"\n****ERROR! Right hand side vector name must be entered\n";
				return 0;
			}
			option->rhs_name = argv[i];
		}
		else if (arg == "of" || arg == "outfile")
		{
			if (!argv[++i])
			{
				std::cout<<"\n****ERROR! Output filename must be entered\n";
				return 0;
			}
			option->out_file = argv[i];
		}
		else
		{
			std::cout<<"\n\n**ERROR! Unknown command line parameter: "<<arg<<"\nIgnoring this parameter\n";
		}
	}
	return 1;
}

/* Purpose: This function reads in a matrix stored in matrix market format
   Arguments:
   A: Matrix to be populated
   filename: name of file containing matrix in matrix market format */
void SLIP_read_matrix(SLIP_mat* A, std::string filename)
{
	SLIP_trip* B = new SLIP_trip;	            // Read in triplet form first
	std::string line;
	std::ifstream in_file;
	in_file.open(filename);
	getline(in_file,line);
	std::istringstream line_ss(line);		    // Read in size of matrix and number of nonzeros
	line_ss >> B->m;
	line_ss >> B->n;
	line_ss >> B->nz;
	B->i = new int [B->nz];				        // Allocate Memory for the triplet i, j and x vectors
	B->j = new int [B->nz];
	B->x = SLIP_initialize_mpz_array(B->nz);	// Initialize the input mpz vector

	for (int i = 0; i < B->nz; i++)			    // Read in the values from file
	{
		getline(in_file,line);
		std::istringstream line_ss(line);
		line_ss >> B->i[i];
		line_ss >> B->j[i];
		line_ss >> B->x[i];
		B->i[i] = B->i[i] - 1;			        // Conversion from 1 based to 0 based
		B->j[i] = B->j[i] - 1;
	}
	in_file.close();
	SLIP_trip_to_ccf(B,A);				        // Convert from triplet form to ccf
	delete B;
}

/* Purpose: This function reads in a RHS vector stored as a dense vector
   This function assumes b is n*1
   Arguments:
   b: mpz array which will store the right hand side vector 
   n: size of b
   filename: file containing b stored as a dense vector */
void SLIP_read_rhs(mpz_t** b, int n, std::string filename)
{
	std::string line;
	std::ifstream in_file;
	in_file.open(filename);
	for (int i = 0; i < n; i++)
	{
		getline(in_file, line);
		std::istringstream line_ss(line);
		line_ss >> b[i][0];
	}
	in_file.close();
}
#endif
