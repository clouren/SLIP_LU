function option = SLIP_get_options
%SLIP_GET_OPTIONS: set defaults for SLIP_LU factorization.
%
% Usage:
%
%   option = SLIP_get_options;
%
% Each element of the option struct has the following possible values:
% option.order: Column ordering used.
%        0: None: not recommended unless you know the matrix already has a
%           good column ordering.
%        1: COLAMD (the default ordering)
%        2: AMD
%
% option.pivot: Row pivoting scheme used.
%        0: Smallest pivot
%        1: Diagonal pivoting
%        2: First nonzero per column chosen as pivot
%        3: Diagonal pivoting with tolerance for smallest pivot (default)
%        4: Diagonal pivoting with tolerance for largest pivot
%        5: Largest pivot
%
% option.tol: tolerance (0,1) which is used if some threshold pivoting is used.
%       default is 0.1.
%
% If any fields are not present in the options struct, the defaults are used
% for that parameter.
%
% See also SLIP_LU, SLIP_install, SLIP_test.

% SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

option.order = 1 ;               % type of column ordering used
option.pivot = 3 ;               % type of row pivoting scheme
option.tol = 0.1 ;               % tolerance for row pivoting scheme

% TODO: if struct fields empty, use the defaults

% TODO: delete this function.

% TODO: use strings, not integers for:
% option.order = 'none', 'colamd', or 'amd'
% option.pivot = 'smallest' , etc

% TODO: delete these options:
option.A_is_integral = false ;   % true if A known to be integral
option.b_is_integral = false ;   % true if b known to be integral

