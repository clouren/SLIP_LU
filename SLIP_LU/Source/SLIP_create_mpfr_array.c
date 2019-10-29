# include "SLIP_LU_internal.h"

/* Purpose: This function creates a MPFR array of desired precision*/
mpfr_t* SLIP_create_mpfr_array
(
    int32_t n,     // size of the array
    SLIP_options *option// command options containing the prec for mpfr
)
{
    // Check input
    if (n <= 0) {return NULL;}
    mpfr_t* x = (mpfr_t*) SLIP_calloc(n, SIZE_MPFR);
    if (!x) {return NULL;}
    for (int32_t i = 0; i < n; i++)        
    {
        if (slip_mpfr_init2(x[i], option->prec) != SLIP_OK)
        {
	    SLIP_MPFR_SET_NULL(x[i]);
            SLIP_delete_mpfr_array(&x, n);
            return NULL;
        }
    }
    return x;
}
