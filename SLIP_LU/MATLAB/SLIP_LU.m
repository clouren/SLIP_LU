%------------------------------------------------------------------------------
% SLIP_LU/MATLAB/SLIP_LU.m: Interface to SLIP LU within MATLAB
%------------------------------------------------------------------------------

% SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

%------------------------------------------------------------------------------

function varargout = SLIP_LU(A,b,option)
% Purpose: Solve exactly the sparse linear system Ax = b where A and b are stored as
% doubles. A must be stored as a sparse matrix. b must be stored as a dense
% set of right hand side vectors. b can be either 1 or multiple vector(s)
%
% ****WARNING****: If A is very large or dense, this function may crash.
%  
% USAGE:
% x = SLIP_LU(A,b) returns the solution to Ax=b using default settings
%
% x = SLIP_LU(A,b,options) returns the solution to Ax=b with user defined
% settings. The options settings can be obtained from option =
% SLIP_get_options then changed from there
%
% [L U P Q x] = SLIP_LU(A,b) returns the solution to the system and lower
% and upper triangular factors L and U such that L*U = P*A*Q using default
% parameters.
%
% [L U P Q x] = SLIP_LU(A,b,options) returns the solution to the system and lower
% and upper triangular factors L and U such that L*U = P*A*Q using user specified
% parameters.
%
% [L U P Q] = SLIP_LU(A) returns lower and upper triangular factors L and 
% U such that L*U = P*A*Q using default parameters.
%
% [L U P Q] = SLIP_LU(A, options) returns lower and upper triangular 
% factors L and U such that L*U = P*A*Q using user specified parameters.

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
elseif (nargout == 5) % x = A\b, L U = PAQ
   option.pivot = 5;
   [varargout{1} varargout{2} varargout{3} varargout{4} ...
       varargout{5}] = SLIP_mex_soln2(A,b,option);
   % Get P and Q
   [m n] = size(varargout{1});
   varargout{3} = varargout{3}+1; varargout{4} = varargout{4}+1;
   A2 = speye(n,n);
   % Return P and Q as permutations of I
   varargout{3} = A2(varargout{3},:);
   varargout{4} = A2(:,varargout{4});
elseif (nargout == 4) % LU = PAQ
    if (nargin > 2)
        error('Incorrect number of input arguments for [L U P Q]. Please type help SLIP_LU')
    end
    if (exist('b') ~= 0)
            if (isstruct(b) == 1)
                option = b;
            end
    end
    option.pivot = 5;
    [varargout{1} varargout{2} varargout{3} varargout{4}] ...
      = SLIP_mex_soln3(A, option);
   % Get P and Q
   [m n] = size(varargout{1});
   varargout{3} = varargout{3}+1; varargout{4} = varargout{4}+1;
   A2 = speye(n,n);
   varargout{3} = A2(varargout{3},:);
   varargout{4} = A2(:,varargout{4});
else
fprintf('Incorrect number of output arguments. Please type help SLIP_LU\n')
end
end
