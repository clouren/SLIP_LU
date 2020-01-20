function maxerr = slip_run_all_tests (A, b)
%SLIP_RUN_ALL_TESTS run a set of sets of SLIP_LU
%
% maxerr = slip_run_all_tests (A, b)

% Test SLIP_mex_soln.c
x = SLIP_LU(A,b);
x2 = A\b;
err1 = norm(x-x2)/norm(x);

maxerr = err1;
