function maxerr = slip_run_all_tests (A, b)
%SLIP_RUN_ALL_TESTS run a set of sets of SLIP_LU
%
% maxerr = slip_run_all_tests (A, b)

% Test SLIP_mex_soln.c
x = SLIP_LU(A,b);
x2 = A\b;
err1 = norm(x-x2)/norm(x);

% Test SLIP_mex_soln2.c
[L U P Q x] = SLIP_LU(A, b);
x2 = A\b;
err2 = norm(x-x2)/(norm(x));
err2a = normest(L*U-P*A*Q);

% Test SLIP_mex_soln3.c
[L U P Q] = SLIP_LU(A);
err3 = normest(L*U-P*A*Q);

maxerr = max ([err1 err2 err2a err3]) ;
