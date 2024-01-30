/**
 * @file findroot.hpp
 * @brief This includes definitions and implementations of 1D-root-finders.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef MATH_FINDROOT_1D_HPP
#define MATH_FINDROOT_1D_HPP

#include <cmath>
#include <iostream>

namespace Math
{
namespace FindRoot1D
{

/**
 * @brief This finds the root of inLowerBound function using the Brent method.
 * @tparam Func_ type of Funtion inF.
 * @param inF function double(double) to find the root
 * @param inLowerBound lower bound of search section
 * @param inUpperBound upper bound of search section
 * @param inTol allowable limit of error
 * @return double the value of root
 */
template < typename Func_ >
double Brent( Func_ inF, double inLowerBound, double inUpperBound,
              double inTol = 1e-6 )
{
    double lFa = inF( inLowerBound );
    double lFb = inF( inUpperBound );

    if ( lFa * lFb >= 0 )
    {
        std::cerr
            << "Error: Math::FindRoot1D::Brent<Func_>(Func_, double, double)"
            << std::endl
            << "The function must have different signs at LowerBound and "
               "UpperBound."
            << std::endl;
        return std::numeric_limits< double >::quiet_NaN();
    }

    if ( std::abs( lFa ) < std::abs( lFb ) )
    {
        std::swap( inLowerBound, inUpperBound );
        std::swap( lFa, lFb );
    }

    double lC   = inLowerBound;
    double lFc  = lFa;
    bool lMflag = true;
    double lS   = inUpperBound;
    double lFs  = lFb;

    while ( lFs != 0 && std::abs( inUpperBound - inLowerBound ) > inTol )
    {
        if ( lFa != lFc && lFb != lFc )
        {
            // Inverse quadratic interpolation
            lS = inLowerBound * lFb * lFc / ( ( lFa - lFb ) * ( lFa - lFc ) ) +
                 inUpperBound * lFa * lFc / ( ( lFb - lFa ) * ( lFb - lFc ) ) +
                 lC * lFa * lFb / ( ( lFc - lFa ) * ( lFc - lFb ) );
        }
        else
        {
            // Secant method
            lS = inUpperBound -
                 lFb * ( inUpperBound - inLowerBound ) / ( lFb - lFa );
        }

        double lTmp1 = ( 3 * inLowerBound + inUpperBound ) / 4;
        double lTmp2 = ( inUpperBound - inLowerBound ) / 2;

        if ( ( !( ( lS > lTmp1 && lS < inUpperBound ) ||
                  ( lS < lTmp1 && lS > inUpperBound ) ) ) ||
             ( lMflag &&
               std::abs( lS - inUpperBound ) >= std::abs( lTmp2 ) / 2 ) ||
             ( !lMflag && std::abs( lS - inUpperBound ) >=
                              std::abs( lC - inUpperBound ) / 2 ) ||
             ( lMflag && std::abs( lTmp2 ) < inTol ) ||
             ( !lMflag && std::abs( lC - inUpperBound ) < inTol ) )
        {
            lS     = ( inLowerBound + inUpperBound ) / 2;
            lMflag = true;
        }
        else { lMflag = false; }

        lFs = inF( lS );
        lC  = inUpperBound;
        lFc = lFb;

        if ( lFa * lFs < 0 )
        {
            inUpperBound = lS;
            lFb          = lFs;
        }
        else
        {
            inLowerBound = lS;
            lFa          = lFs;
        }

        if ( std::abs( lFa ) < std::abs( lFb ) )
        {
            std::swap( inLowerBound, inUpperBound );
            std::swap( lFa, lFb );
        }
    }
    return inUpperBound;
}

}  // namespace FindRoot1D
}  // namespace Math

#endif