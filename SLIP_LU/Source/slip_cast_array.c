//------------------------------------------------------------------------------
// SLIP_LU/slip_cast_array: scale and typecast an array
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothx A. Davis, Teyas A&M Universitx.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

// Scales and typecasts an input array X, into the output array Y.

// X: an array of type xtype, of size n.
// Y: an array of type ytype, of size n.

// Y [0:n-1] = scale * X [0:n-1] ;

#include "SLIP_LU_internal.h"
#pragma GCC diagnostic ignored "-Wunused-variable"

SLIP_info slip_cast_array
(
    void *Y,                // output array, of size n
    SLIP_type ytype,        // type of Y
    void *X,                // input array, of size n
    SLIP_type xtype,        // type of X
    int32_t n,              // size of Y and X
    mpq_t scale,            // scale factor applied if X is mpz_t
    SLIP_options *option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SLIP_info info ;
    if (n <= 0)
    {
        // nothing to do
        return (SLIP_OK) ;
    }
    if (Y == NULL || X == NULL || option == NULL)
    {
        return (SLIP_INCORRECT_INPUT) ;
    }

    //--------------------------------------------------------------------------
    // Y [0:n-1] = (ytype) X [0:n-1]
    //--------------------------------------------------------------------------

    switch (ytype)
    {

        //----------------------------------------------------------------------
        // output array Y is mpz_t
        //----------------------------------------------------------------------

        case SLIP_MPZ:
        {
            mpz_t *y = (mpz_t *) Y ;
            switch (xtype)
            {

                case SLIP_MPZ: // mpz_t to mpz_t
                {
                    mpz_t *x = (mpz_t *) X ;
                    for (int32_t k = 0 ; k < n ; k++)
                    {
                        SLIP_CHECK (SLIP_mpz_set (y [k], x[k])) ;
                    }
                }
                break ;

                case SLIP_MPQ: // mpq_t to mpz_t
                {
                    mpq_t *x = (mpq_t *) X ;
                    SLIP_CHECK (slip_expand_mpq_array (Y, X, scale, n)) ;
                }
                break ;

                case SLIP_MPFR: // mpfr_t to mpz_t
                {
                    mpfr_t *x = (mpfr_t *) X ;
                    SLIP_CHECK (slip_expand_mpfr_array (Y, X, scale, n,
                        option)) ;
                }
                break ;

                case SLIP_INT32: // int32_t to mpz_t
                {
                    int32_t *x = (int32_t *) X ;
                    for (int32_t k = 0 ; k < n ; k++)
                    {
                        SLIP_CHECK (SLIP_mpz_set_si (y [k], x [k])) ;
                    }
                    SLIP_CHECK (SLIP_mpq_set_ui (scale, 1, 1)) ;
                }
                break ;

                case SLIP_FP64: // double to mpz_t
                {
                    double *x = (double *) X ;
                    SLIP_CHECK (slip_expand_double_array (y, x, scale, n,
                        option)) ;
                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
            }
        }
        break ;

        //----------------------------------------------------------------------
        // output array Y is mpq_t
        //----------------------------------------------------------------------

        case SLIP_MPQ:
        {
            mpq_t *y = (mpq_t *) Y ;
            switch (xtype)
            {

                case SLIP_MPZ: // mpz_t to mpq_t
                {
                    mpz_t *x = (mpz_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_MPQ: // mpq_t to mpq_t
                {
                    mpq_t *x = (mpq_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_MPFR: // mpfr_t to mpq_t
                {
                    mpfr_t *x = (mpfr_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_INT32: // int32 to mpq_t
                {
                    int32_t *x = (int32_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_FP64: // double to mpq_t
                {
                    double *x = (double *) X ;
                    // TODO
                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
            }
        }
        break ;

        //----------------------------------------------------------------------
        // output array Y is mpfr_t
        //----------------------------------------------------------------------

        case SLIP_MPFR:
        {
            mpfr_t *y = (mpfr_t *) Y ;
            switch (xtype)
            {

                case SLIP_MPZ: // mpz_t to mpfr_t
                {
                    mpz_t *x = (mpz_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_MPQ: // mpq_t to mpfr_t
                {
                    mpq_t *x = (mpq_t *) X ;
                    for (int32_t k = 0 ; k < n ; k++)
                    {
                        SLIP_CHECK (SLIP_mpfr_set_q (y [k], x [k],
                            option->SLIP_MPFR_ROUND)) ;
                    }
                }
                break ;

                case SLIP_MPFR: // mpfr_t to mpfr_t
                {
                    mpfr_t *x = (mpfr_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_INT32: // int32 to mpfr_t
                {
                    int32_t *x = (int32_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_FP64:  // double to mpfr_t
                {
                    double *x = (double *) X ;
                    // TODO
                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
            }
        }
        break ;

        //----------------------------------------------------------------------
        // output array Y is int32_t
        //----------------------------------------------------------------------

        case SLIP_INT32:
        {
            int32_t *y = (int32_t *) Y ;
            switch (xtype)
            {

                case SLIP_MPZ: // mpz_t to int32_t
                {
                    mpz_t *x = (mpz_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_MPQ: // mpq_t to int32_t
                {
                    mpq_t *x = (mpq_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_MPFR: // mpfr_t to int32_t
                {
                    mpfr_t *x = (mpfr_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_INT32: // int32_t to int32_t
                {
                    memcpy (Y, X, n * sizeof (int32_t)) ;
                }
                break ;

                case SLIP_FP64: // double to int32_t
                {
                    double *x = (double *) X ;
                    for (int32_t k = 0 ; k < n ; k++)
                    {
                        if (isnan (x [k]))
                        {
                            y [k] = 0 ;
                        }
                        else if (x [k] > INT32_MAX)
                        {
                            y [k] = INT32_MAX ;
                        }
                        else if (x [k] < INT32_MIN)
                        {
                            y [k] = INT32_MIN ;
                        }
                        else
                        {
                            y [k] = (int32_t) (x [k]) ;
                        }
                    }
                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
            }
        }
        break ;

        //----------------------------------------------------------------------
        // output array Y is double
        //----------------------------------------------------------------------

        case SLIP_FP64:
        {
            double *y = (double *) Y ;
            switch (xtype)
            {

                case SLIP_MPZ: // mpz_t to double
                {
                    mpz_t *x = (mpz_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_MPQ: // mpq_t to double
                {
                    mpq_t *x = (mpq_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_MPFR: // mpfr_t to double
                {
                    mpfr_t *x = (mpfr_t *) X ;
                    // TODO
                }
                break ;

                case SLIP_INT32: // int32_t to double
                {
                    int32_t *x = (int32_t *) X ;
                    for (int32_t k = 0 ; k < n ; k++)
                    {
                        y [k] = (double) (x [k]) ;
                    }
                }
                break ;

                case SLIP_FP64: // double to double
                {
                    memcpy (Y, X, n * sizeof (double)) ;
                }
                break ;

                default: return (SLIP_INCORRECT_INPUT) ;
            }
        }
            break ;

        default: return (SLIP_INCORRECT_INPUT) ;
    }

    return (SLIP_OK) ;
}
