function x = SLIP_LU(A,b,option)
%SLIP_LU: solve Ax=b via sparse left-looking integer-preserving LU factorization
%
% Purpose: Solve exactly the sparse linear system Ax = b where A and b are
% stored as doubles. A must be stored as a sparse matrix. b must be stored as a
% dense set of right hand side vectors. b can be either 1 or multiple vector(s)
%
% ****WARNING****: If A is very large or dense, this function may crash.
% TODO: how?  Will it segfault and MATLAB dies?  Or just report out-of-memory?
% If the latter, then delete this statement (it is true for all MATLAB
% functions).  If the former ... why?  It shouldn't do that.
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

% Check for integer overflow. If the max value of A or b exceeds this value
% the internal routines can not expect integer input.

% TODO what is this magic number 2000000000000 ??
% It's 2e12, which is a lot bigger than 2^31.

if (max(max(abs(A))) > 2000000000000)
    option.A_is_integral = false ;
end
if (max(abs(A)) > 2000000000000)
    option.b_is_integral = false ;
end

% Check if the input matrix is stored as sparse. If not, SLIP LU expects
% sparse input, so convert to sparse.
if (~issparse (A))
    A = sparse (A) ;
end

% If the user indicates that the input is integral, check if it is actually
% integral.  TODO: do this inside the mexFunction.

if (option.A_is_integral)
    A2 = floor(A);
    % TODO: so if it is < 1e-12, then the matrix can be treated as 'integeral'?
    % Doesn't that introduce some error?  Why not an exact check?
    if (normest(A2-A) > 1e-12)
        option.A_is_integral = false ;
    end
    clear A2;
end

if (option.b_is_integral)
    b2 = floor(b);
    if (normest(b2-b) > 1e-12)
        option.b_is_integral = false ;
    end
    clear b2;
end

% Preprocessing complete. Now use SLIP LU to solve Ax=b.
x = SLIP_mex_soln (A, b, option) ;

