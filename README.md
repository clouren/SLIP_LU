SLIP_LU is a software package for solving a sparse systems of linear
equations exactly using the SLIP LU factorization. 

Files and folders in this distribution:

    README.md   this file
    SLIP_LU     the primary SLIP_LU library, demos, and test programs
    Makefile    compiles SLIP_LU and its dependencies

Dependencies (all part of SuiteSparse):

    AMD         
    COLAMD
    SuiteSparse_config

Default instalation locations:

    include
    lib
    share

To compile SLIP_LU and its dependencies, just type "make" in this folder.
This will also run a few short demos in `SLIP_LU/Demo`.
To install the package system-wide, copy the `lib/*` to /usr/local/lib,
and copy `include/*` to /usr/local/include.

