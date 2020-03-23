function x = SLIP_LU(A,b,option)
%SLIP_LU: solve Ax=b via sparse left-looking integer-preserving LU factorization
%
% Purpose: Solve exactly the sparse linear system Ax = b where A and b are
% stored as doubles. A must be stored as a sparse matrix. b must be stored as a
% dense set of right hand side vectors. b can be either 1 or multiple vector(s)
%
% USAGE:
%
% x = SLIP_LU(A,b) returns the solution to Ax=b using default settings
%
% x = SLIP_LU(A,b,options) returns the solution to Ax=b with user defined
% settings. The options settings can be obtained from option = SLIP_get_options
% then changed from there.
%
% See also SLIP_install, SLIP_get_options, SLIP_test.

% SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

if (nargin < 3)
    option = SLIP_get_options ;   % Set defaults
end

% Check if the input matrix is stored as sparse. If not, SLIP LU expects
% sparse input, so convert to sparse.
if (~issparse (A))
    A = sparse (A) ;
end

% Preprocessing complete. Now use SLIP LU to solve Ax=b.
x = SLIP_mex_soln (A, b, option) ;

