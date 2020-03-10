function varargout = SLIP_LU(A,b,option)
%SLIP_LU: solve Ax=b via sparse left-looking integer-preserving LU factorization
%
% Purpose: Solve exactly the sparse linear system Ax = b where A and b are
% stored as doubles. A must be stored as a sparse matrix. b must be stored as a
% dense set of right hand side vectors. b can be either 1 or multiple vector(s)
%
% ****WARNING****: If A is very large or dense, this function may crash.
%
% USAGE:
%
% x = SLIP_LU(A,b) returns the solution to Ax=b using default settings
%
% x = SLIP_LU(A,b,options) returns the solution to Ax=b with user defined
% settings. The options settings can be obtained from option = SLIP_get_options
% then changed from there.

% SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

if exist('option') == 0          % Did the user pass in options?
    option = SLIP_get_options;   % Set defaults
end

% Check for int overflow. If the max value of A or b exceeds this value
% the internal routines can not expect int input.
if (max(max(abs(A))) > 2000000000000)
    option.int = 0;
end
if (max(abs(A)) > 2000000000000)
    option.intb = 0;
end

% Check if the input matrix is stored as sparse. If not, SLIP LU expects
% sparse input, so convert to sparse.
if (issparse(A) == 0)
    A = sparse(A);
end

% If the user indicates that the input is integral, check if it is actually integral
if (option.int > 0)
    A2 = floor(A);
    if (normest(A2-A) > 1e-12)
        option.int = 0;
    end
    clear A2;
end

if (option.intb > 0)
    b2 = floor(b);
    if (normest(b2-b) > 1e-12)
        option.intb = 0;
    end
    clear b2;
end

% Preprocessing complete. Now use SLIP LU
if (nargout == 1) % x = A\b
    varargout{1} = SLIP_mex_soln(A,b,option);
else
    fprintf('Incorrect number of output arguments. Please type help SLIP_LU\n')
end

