/**
 * @file special_functions.hpp
 * @brief This defines special functions.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef MATH_SPECIAL_FUNCTIONS_HPP
#define MATH_SPECIAL_FUNCTIONS_HPP

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

}  // namespace SpecialFunctions
}  // namespace Math

#endif