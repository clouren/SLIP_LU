%%SLIP_DEMO a demo of SLIP_backslash
%
% See also SLIP_backslash, SLIP_install, SLIP_test.
%
% SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
% SLIP_LU/License for the license.

%% SLIP_backslash vs MATLAB backslash: first example
% In this first example, x = SLIP_backslash (A,b) returns an approximate
% solution, not because it was computed incorrectly in SLIP_backslash.  It
% is computed exactly as a rational result in SLIP_backslash with arbitrary
% precision, but then converted to double precision on output.

load west0479
A = west0479 ;
n = size (A, 1) ;
xtrue = rand (n,1) ;
b = A*xtrue ;
x = SLIP_backslash (A, b) ;
% error is nonzero: x is computed exactly in rational arbitrary-precision,
% but then lost precision when returned to MATLAB:
err_slip = norm (x-xtrue)
x = A\b ;
% error is nonzero: MATLAB x=A\b experiences floating-point error
% throughout its computations:
err_matlab = norm (x-xtrue)

%% SLIP_backslash: exact, vs MATLAB backslash: approximate
% In this example, x = SLIP_backslash (A,b) is returned exactly in the
% MATLAB vector x, because x contains only integers representable exactly
% in double precision.  x = A\b results in floating-point roundoff error.

amax = max (abs (A), [ ], 'all') ;
A = floor (2^20 * (A / amax)) + n * speye (n) ;
xtrue = floor (64 * xtrue) ;
b = A*xtrue ;
x = SLIP_backslash (A, b) ;
% error will be exactly zero:
err_slip = norm (x-xtrue)
x = A\b ;
% error will be small but nonzero:
err_matlab = norm (x-xtrue)

%% SLIP_backslash on difficult problems

% TODO add some sample solutions from some nasty
% matrices from the MATLAB gallery function.

