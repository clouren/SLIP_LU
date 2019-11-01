% This function installs the SLIP LU matlab routines. It allows the use of
% m files SLIP_LU.m and SLIP_det.m
% Please run this command prior to attempting to use any SLIP LU routines
% Usage: SLIP_install
% Required Libraries: GMP, MPFR, AMD, COLAMD

% Find all source files and add them to the src string
src = '';
path = './Source/';
files = dir('./Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end
path = '../Source/';
files = dir('../Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end

% Compiler flags
flags = 'CFLAGS=''-std=c99 -fPIC''';

% External libraries
libs = '-lgmp -lmpfr -lamd -lcolamd';

% Path to headers
includes = '-ISource/ -I../Source/ -I../Include/ ';

% Generate the mex commands here
% having -R2018a here for function mxGetDoubles
m1 = ['mex -R2018a ', includes, ' SLIP_mex_soln.c ' , src, ' ', flags, ' ', libs];
m2 = ['mex -R2018a ', includes, ' SLIP_mex_soln2.c ', src, ' ', flags, ' ', libs];
m3 = ['mex -R2018a ', includes, ' SLIP_mex_soln3.c ', src, ' ', flags, ' ', libs];

% Now, we evaluate each one
eval(m1);
eval(m2);
eval(m3);

fprintf('\nMex files installed, now we test\n')

% Efficient testing
SLIP_test;
