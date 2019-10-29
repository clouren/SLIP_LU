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
% [x d] = SLIP_LU(A,b) returns the solution to the system and the absolute
% value of the determinant of A
%
% [x d] = SLIP_LU(A,b,options) returns the solution to the system and the
% absolute value of the determinant of A using user defined options.
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
    option = SLIP_get_options;  % Set defaults
end

if (max(max(abs(A))) > 2000000000000) % Integer overflow
    option.int = 0;
end
if (max(abs(A)) > 2000000000000)
    option.intb = 0;
end
if (issparse(A) == 0)
    A = sparse(A);
end

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

if (nargout == 1) % x = A\b
    varargout{1} = SLIP_mex_soln(A,b,option);
elseif (nargout == 2) % x = A\b, d = det(A)
   [varargout{1} varargout{2}] = SLIP_mex_soln(A,b,option);
elseif (nargout == 5) % x = A\b, L U = PAQ
   option.pivot = 5;
   [varargout{1} varargout{2} varargout{3} varargout{4} ...
       varargout{5}] = SLIP_mex_soln2(A,b,option);
   % Get P and Q
   [m n] = size(varargout{1});
   varargout{3} = varargout{3}+1; varargout{4} = varargout{4}+1;
   p2 = zeros(n,1);
   A2 = speye(n,n);
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
   p2 = zeros(n,1);
   A2 = speye(n,n);
   varargout{3} = A2(varargout{3},:);
   varargout{4} = A2(:,varargout{4});
else
fprintf('Incorrect number of output arguments. Please type help SLIP_LU\n')
end
end
