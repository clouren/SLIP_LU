//------------------------------------------------------------------------------
// SLIP_LU/SLIP_get_double_soln: convert mpq solution to double
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Convert the output mpq_t** solution vector obtained from SLIP_Solve
 * and SLIP_Permute_x from mpq_t** to double.  x_doub has to be initialized
 * before passed in.
 *
 * Input/output arguments:
 *
 * x_doub:  double** containing the users final solution in double precision.
 *          uninitialized on input, contains the double version of the solution
 *          on output.
 *
 * x_mpq:   mpq_t** containing the exact solution of the linear system.
 *          Unmodified on output and input.
 *
 * n:       Number of columns in the input matrix = number of rows in x
 *
 * numRHS:  Number of RHS vectors = number of columns in x
 */

#include "SLIP_LU_internal.h"

SLIP_info SLIP_get_double_soln
(
    double **x_doub,      // double soln of size n*numRHS to Ax = b
    mpq_t  **x_mpq,       // mpq solution to Ax = b. x is of size n*numRHS
    int32_t n,            // Dimension of A, number of rows of x
    int32_t numRHS        // Number of right hand side vectors
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    // TODO: use SLIP_matrix_copy (&X_double, SLIP_DENSE, SLIP_DOUBLE, X_mpq...)
    // and delete this function.

    SLIP_info info ;
    if (x_doub  == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------

    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t j = 0; j < numRHS; j++)
        {
            SLIP_CHECK (SLIP_mpq_get_d (&(x_doub[i][j]), x_mpq[i][j])) ;
        }
    }
    return SLIP_OK;
}

