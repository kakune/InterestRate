/**
 * @file special_functions.hpp
 * @brief This includes definitions and implementations of special functions.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef MATH_SPECIAL_FUNCTIONS_HPP
#define MATH_SPECIAL_FUNCTIONS_HPP

#include <cmath>

namespace Math
{

namespace SpecialFunctions
{

/**
 * @brief This calculates CDF of normal distribution.
 * @param inX argument of CDF
 * @return double normal CDF N(inX)
 */
double normalCDF( double inX );

/**
 * @brief This calculates PDF of normal distribution.
 * @param inX argument of PDF
 * @return double normal PDF \phi(inX)
 */
double normalPDF( double inX );

double normalCDF( double inX ) { return 0.5 * std::erfc( -inX * M_SQRT1_2 ); }
double normalPDF( double inX )
{
    return ( 1.0 / std::sqrt( 2.0 * M_PI ) ) * std::exp( -0.5 * inX * inX );
}

}  // namespace SpecialFunctions

}  // namespace Math

#endif