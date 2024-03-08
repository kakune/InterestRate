/**
 * @file integral_1d.hpp
 * @brief This defines the 1d-integral functions
 * @author kakune
 * @date 3/8/2024
 */

#ifndef MATH_INTEGRAL_1D
#define MATH_INTEGRAL_1D

#include <vector>

namespace Math
{
namespace Integral
{

std::vector<double> eachTrapezoidal( const std::vector<double>& inXs,
                                     const std::vector<double>& inYs );

}
}  // namespace Math

#endif