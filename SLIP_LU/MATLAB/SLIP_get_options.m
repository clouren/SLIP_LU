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
%
% option.pivot: Pivoting scheme used.
%        0: Smallest pivot
%        1: Diagonal pivoting
%        2: First nonzero per column chosen as pivot
%        3: Diagonal pivoting with tolerance for smallest pivot, Default
%        4: Diagonal pivoting with tolerance for largest pivot
%        5: Largest pivot
% 
% option.int: Set equal to 1 if input mat is already integral
% option.intb: set equal to 1 if input RHS vector(s) are already integral
% option.tol: tolerance (0,1) which is used if some threshold pivoting is used

x.column = 1;
x.pivot = 3;
x.int = 0;
x.intb = 0;
x.tol = 0.1;
end
