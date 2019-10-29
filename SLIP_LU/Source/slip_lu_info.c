//------------------------------------------------------------------------------
// SLIP_LU/slip_lu_info: print the version and date of SLIP_LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_internal.h"

void slip_lu_info()
{
    printf("\n****Software Information****");
    printf("\n\nYou are using SLIP_LU Version: %s",SLIP_LU_VERSION);
    printf("\nThis code accompanies the paper:\n\n\t%s\n",SLIP_PAPER);
    printf("\nThis code is copyright by %s\n\n",SLIP_AUTHOR);
}
