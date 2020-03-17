function option = SLIP_get_options
%SLIP_GET_OPTIONS: set defaults for SLIP_LU factorization.
%
% Usage:
%
%   option = SLIP_get_options ;
%   option.tol = 1e-6 ;
%   x = SLIP_LU (A, b, option) ;
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
% If a field is not present in the options struct, or if the options
% struct is not passed to SLIP_LU, the defaults are used.
%
% See also SLIP_LU, SLIP_install, SLIP_test.

% SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

option.order = 1 ;               % type of column ordering used
option.pivot = 3 ;               % type of row pivoting scheme
option.tol = 0.1 ;               % tolerance for row pivoting scheme

% TODO: delete this function; not needed

% TODO: use strings, not integers for:
% option.order = 'none', 'colamd', or 'amd'
% option.pivot = 'smallest' , etc

