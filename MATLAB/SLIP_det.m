function d = SLIP_det(A, option)
% Purpose: Compute the determinant of A. A must be stored as a sparse matrix.
%
% ****WARNING****: If A is very large or dense, this function may crash or have severely
% long run time.
%
% USAGE:
% d = SLIP_det(A) returns the determinant of A using default settings
%
% d = SLIP_det(A,option) returns the determinant of A using 
% user defined settings
if exist('option') ==0          % Did the user pass in options?
    option = SLIP_get_options;  % Set defaults
end

if (max(max(abs(A))) > 2000000000000) % Integer overflow
    option.int = 0;
end
if (issparse(A) == 0) % Is A sparse?
    A = sparse(A);
end

if (option.int > 0) % Does the user think A is integral?
    A2 = floor(A);
    if (normest(A2-A) > 1e-12) % Is INT(A) == A?
        option.int = 0;
    end
    clear A2;
end

if (nargout == 1) % x = A\b
    d = SLIP_mex_soln4(A,option);
else
fprintf('Usage: d = SLIP_det(A,option)\n')
end
end
