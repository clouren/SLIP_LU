function x = SLIP_get_options()
% This function defines various parameters for the SLIP LU factorization. It sets
% each parameter to its default value.
% 
% USAGE: option = SLIP_get_options;
%
% Each element of the option struct has the following possible values:
% option.column: Column ordering used. 
%        0: None: Not recommended for sparse matrices
%        1: COLAMD: Default
%        2: AMD
%        3: UMFPACK P and Q
%
% option.pivot: Pivoting scheme used. Note that if UMFPACK is selected for
% column, UMFPACK's pivot order is automatically used.
%        0: Smallest pivot: Default and recommended
%        1: Diagonal pivoting
%        2: First nonzero per column chosen as pivot
%        3: Diagonal pivoting with tolerance for smallest pivot
%        4: Diagonal pivoting with tolerance for largest pivot
%        5: Largest pivot
% 
% option.print_level: level of details to be printed
%        0: print nothing
%        1: just errors and warnings: Default
%        2: terse, with basic stats from COLAMD/AMD and SLIP and solution
%
% option.check: Set equal to 1 if you want solution checked 
% option.int: Set equal to 1 if input mat is already integral
% option.intb: set equal to 1 if input RHS vector(s) are already integral
% option.tol: tolerance (0,1) which is used if some threshold pivoting is used

x.column = 1;
x.pivot = 0;
x.print_level = 0;
x.check = 0;
x.int = 0;
x.intb = 0;
x.tol = 0.1;
end