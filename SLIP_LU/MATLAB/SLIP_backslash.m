function x = SLIP_backslash (A,b,option)
%SLIP_BACKSLASH: solve Ax=b via sparse left-looking integer-preserving LU
% SLIP_backslash: computes the exact solution to the sparse linear system Ax =
% b where A and b are stored as doubles. A must be stored as a sparse matrix. b
% must be stored as a dense set of right hand side vectors. b can be either 1
% or multiple vector(s).  The result x is computed exactly, represented in
% arbitrary-precision rational values, and then returned to MATLAB as a
% floating-poing double result.  This final conversion means that x may no
% longer exactly solve A*x=b, unless this final conversion is able to be
% done without modification.
%
% Usage:
%
% x = SLIP_backslas (A,b) returns the solution to Ax=b using default settings
%
% x = SLIP_backslas (A,b,options) returns the solution to Ax=b with user defined
% settings in an options struct.  Entries not present are treated as defaults.
%
%   option.order: Column ordering used.
%       'none': no column ordering; factorize the matrix A as-is
%       'colamd': COLAMD (the default ordering)
%       'amd': AMD
%
%   option.pivot: Row pivoting scheme used.
%       'smallest': Smallest pivot
%       'diagonal': Diagonal pivoting
%       'first': First nonzero per column chosen as pivot
%       'tol smallest': Diagonal pivoting with tol for smallest pivot (default)
%       'tol largest': Diagonal pivoting with tol for largest pivot
%       'largest': Largest pivot
%
%   option.tol: tolerance (0,1] for 'tol smallest' or 'tol largest' pivoting.
%       default is 0.1.
%
% Example:
%
%   % In this first example, x = SLIP_backslash (A,b) returns an approximate
%   % solution, not because it was computed incorrectly in SLIP_backslash.  It
%   % is computed exactly as a rational result in SLIP_backslash with arbitrary
%   % precision, but then converted to double precision on output.
%
%   load west0479
%   A = west0479 ;
%   n = size (A, 1) ;
%   xtrue = rand (n,1) ;
%   b = A*xtrue ;
%   x = SLIP_backslash (A, b) ;
%   err = norm (x-xtrue)
%   x = A\b ;
%   err = norm (x-xtrue)
%
%   % In this example, x = SLIP_backslash (A,b) is returned exactly in the
%   % MATLAB vector x, because x contains only integers representable exactly
%   % in double precision.  x = A\b results in floating-point roundoff error.
%
%   amax = max (abs (A), [ ], 'all') ;
%   A = floor (2^20 * (A / amax)) + n * speye (n) ;
%   xtrue = floor (64 * xtrue) ;
%   b = A*xtrue ;
%   x = SLIP_backslash (A, b) ;
%   % error and residual will be exactly zero:
%   err = norm (x-xtrue)
%   resid = norm (A*x-b)
%   x = A\b ;
%   % error and residual will be nonzero:
%   err = norm (x-xtrue)
%   resid = norm (A*x-b)
%
% See also SLIP_install, SLIP_test, SLIP_demo.

% SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

if (nargin < 3)
    option = [ ] ;   % use default options
end

% Check if the input matrix is stored as sparse. If not, SLIP LU expects
% sparse input, so convert to sparse.
if (~issparse (A))
    A = sparse (A) ;
end

% Preprocessing complete. Now use SLIP LU to solve Ax=b.
x = SLIP_mex_soln (A, b, option) ;

