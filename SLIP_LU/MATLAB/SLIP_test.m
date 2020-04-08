function SLIP_test
%SLIP_test: run a set of tests for SLIP_LU
%
% Usage:  SLIP_test
%
% See also SLIP_install, SLIP_get_options, SLIP_LU.

maxerr = 0 ;
rng ('default') ;

fprintf ('Testing SLIP_LU: ') ;

% First, check if we can use a real life sparse matrix via ssget
if (exist ('ssget') ~= 0)
    fprintf ('. (please wait) ') ;
    % 159 is a square SPD matrix
    prob = ssget(159);
    A = prob.A;
    [m n] = size(A);
    b = rand(m, 1);
    fprintf ('.') ;
    x = SLIP_LU(A,b);
    x2 = A\b;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;

    % now convert to an integer problem (x will not be integer)
    A = floor (2^20 * A) ;
    b = floor (2^20 * b) ;
    fprintf ('.') ;
    x = SLIP_LU (A, b) ;
    x2 = A\b;
    err = norm(x-x2)/norm(x);
    maxerr = max (maxerr, err) ;
    fprintf ('.') ;
end

for n = [1 10 100]
    for density = [0.001 0.05 0.5 1]
        fprintf ('.') ;
        A = sprand(n,n,density);
        A = A+A';
        % Want a numerically stable A
        if (condest(A) > 1e6)
            A = A + speye (n) ;
        end
        b = rand(n,1);
        x = SLIP_LU(A,b);
        x2 = A\b;
        err = norm(x-x2)/norm(x);
        maxerr = max (maxerr, err) ;

        % now convert to an integer problem (x will not be integer)
        A = floor (2^20 * A) ;
        b = floor (2^20 * b) ;
        x = SLIP_LU(A,b);
        x2 = A\b;
        err = norm(x-x2)/norm(x);
        maxerr = max (maxerr, err) ;

    end
end

fprintf ('\nmaxerr: %g\n', maxerr) ;

if (maxerr < 1e-6)
    fprintf('\nTesting complete, installation successful\n')
else
    error ('SLIP_LU:test', '\nTesting failure!  error too high\n')
end
