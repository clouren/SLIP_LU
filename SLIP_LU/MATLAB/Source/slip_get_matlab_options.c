//------------------------------------------------------------------------------
// SLIP_LU/MATLAB/slip_get_matlab_options: Set factorization options for SLIP LU
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

#include "SLIP_LU_mex.h"

/* Purpose: This function reads in the necessary information from the options
   struct for matlab */
void slip_get_matlab_options
(
    SLIP_options* option,  // Control parameters (must not be NULL)
    const mxArray* input   // The input options from MATLAB interface
)
{
    mxArray* tmp;

    // Get the column ordering
    int64_t order = 1 ;     // default: COLAMD ordering
    tmp = mxGetField(input, 0, "order");
    if (tmp != NULL)
    {
        order = (int64_t) mxGetScalar(tmp);
    }

    // Get the row pivoting scheme
    int64_t piv = 3 ;       // default: diag pivoting with tolerance
    tmp = mxGetField(input, 0, "pivot");
    if (tmp != NULL)
    {
        piv = (int64_t) mxGetScalar(tmp);
    }

    // tolerance for row partial pivoting
    option->tol = 0.1 ;     // default tolerance is 0.1
    if (piv == 3 || piv == 4)
    {
        tmp = mxGetField(input, 0, "tol");
        if (tmp != NULL)
        {
            option->tol = mxGetScalar(tmp);
        }
    }

    //--------------------------------------------------------------------------
    // Verify that the parameters are correct (use defaults if out of range)
    //--------------------------------------------------------------------------

    if (order <= 2 && order >= 0)
    {
        option->order = (SLIP_col_order) order;
    }

    if (piv <= 5 && piv >= 0)
    {
        option->pivot = (SLIP_pivot) piv ;
    }

    if (option->tol > 1 || option->tol <= 0)
    {
        option->tol = 0.1 ;
    }
}

