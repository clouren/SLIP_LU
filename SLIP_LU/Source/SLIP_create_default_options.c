# include "SLIP_LU_internal.h"

/* Purpose: Create and return SLIP_options pointer with default parameters
 * upon successful allocation, which are defined in SLIP_LU_internal.h 
 */

#define CHECK_RESULT(ok)                \
{                                       \
    if (ok != SLIP_OK)                  \
    {                                   \
        SLIP_delete_options(&option);   \
        return NULL;                    \
    }                                   \
}

SLIP_options* SLIP_create_default_options ( void )
{
    SLIP_options* option = SLIP_malloc(sizeof(SLIP_options));
    if (!option) {return NULL;}
    option->check    = SLIP_DEFAULT_CHECK;
    option->pivot    = SLIP_DEFAULT_PIVOT;
    option->order    = SLIP_DEFAULT_ORDER;
    option->print_level= SLIP_DEFAULT_PRINT_LEVEL;
    option->prec     = SLIP_DEFAULT_PRECISION;
    option->tol      = SLIP_DEFAULT_TOL;
    return option;
}
#undef CHECK_RESULT
