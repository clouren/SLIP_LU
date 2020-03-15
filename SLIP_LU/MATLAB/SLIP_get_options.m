function x = SLIP_get_options
%SLIP_GET_OPTIONS: set defaults for SLIP_LU factorization.
%
% Usage:
%
%   option = SLIP_get_options;
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
%
% See also SLIP_LU, SLIP_install, SLIP_test.

% SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

x.column = 1;
x.pivot = 3;
x.int = 0;
x.intb = 0;
x.tol = 0.1;

