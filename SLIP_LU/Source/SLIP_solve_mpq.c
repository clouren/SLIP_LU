# include "SLIP_LU_internal.h"

/* Purpose: This code utilizes the SLIP LU factorization. 
 * Soln is output as mpq_t mat.
 */

# define SLIP_FREE_WORKSPACE        \
    SLIP_delete_sparse(&L);         \
    SLIP_delete_sparse(&U);         \
    SLIP_FREE(pinv);                \
    SLIP_delete_mpz_array(&rhos,n);

SLIP_info SLIP_solve_mpq 
( 
    mpq_t **x_mpq,          // Solution vector stored as an mpq_t array 
    SLIP_sparse *A,         // Compressed column form full precision matrix A 
    SLIP_LU_analysis *S,    // Column ordering  
    SLIP_dense *b,          // Right hand side vectrors 
    SLIP_options *option,   // Control parameters 
    FILE *file              // file to print to, NULL if not used
)
{
    //-------------------------------------------------------------------------
    // Check input
    //-------------------------------------------------------------------------
    if (!x_mpq || !A || !A->p || !A->i || !A->x || 
        !S || !S->q || !b || !b->x || !option)
    {
	return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare memory
    //--------------------------------------------------------------------------
    int32_t *pinv, n = A->n, numRHS = b->n;
    SLIP_info ok, check2 = SLIP_OK ;
    SLIP_sparse* L = SLIP_create_sparse();
    SLIP_sparse* U = SLIP_create_sparse();
    pinv = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
    mpz_t* rhos = SLIP_create_mpz_array(n);
    if (!L || !U || !pinv || !rhos) 
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    /* LU Factorization */
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_LU_factorize(L, U, A, S, rhos, pinv, option));
    
    //--------------------------------------------------------------------------
    // FB Substitution
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_LU_solve(x_mpq, b, rhos, L, U, pinv));
    
    SLIP_CHECK(SLIP_permute_x(x_mpq, n, numRHS, S));
    
    //--------------------------------------------------------------------------
    // Check solution
    //--------------------------------------------------------------------------
    if (option->check)
    {
        SLIP_CHECK(SLIP_check_solution(A, x_mpq, b));
	check2 = ok;
    }
    SLIP_CHECK(SLIP_scale_x(x_mpq, A, b));

    //--------------------------------------------------------------------------
    // Print stats, free memory
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_print_stats_mpq(file, x_mpq, n, numRHS, check2, option));

    SLIP_FREE_WORKSPACE;
    return ok;
}
#undef SLIP_FREE_WORKSPACE
