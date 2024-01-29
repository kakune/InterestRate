/**
 * @file special_functions.cpp
 * @brief This implements special functions.
 * @author kakune
 * @date 1/29/2024
 */

#include <cmath>

namespace Math
{
namespace SpecialFunctions
{

double normalCDF( double inX ) { return 0.5 * std::erfc( -inX * M_SQRT1_2 ); }
double normalPDF( double inX )
{
    return ( 1.0 / std::sqrt( 2.0 * M_PI ) ) * std::exp( -0.5 * inX * inX );
}

}  // namespace SpecialFunctions
}  // namespace Math
