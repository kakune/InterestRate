/**
 * @file integral_1d.cpp
 * @brief This implements the 1d-integral functions
 * @author kakune
 * @date 3/8/2024
 */

#include "math/integral_1d.hpp"

namespace Math
{
namespace Integral
{

std::vector<double> eachTrapezoidal( const std::vector<double>& inXs,
                                     const std::vector<double>& inYs )
{
    std::vector<double> lResult( inXs.size(), 0.0 );
    for ( std::size_t iX = 1; iX < inXs.size(); ++iX )
    {
        lResult.at( iX ) =
            lResult.at( iX - 1 ) + 0.5 * ( inXs.at( iX ) - inXs.at( iX - 1 ) ) *
                                       ( inYs.at( iX ) + inYs.at( iX - 1 ) );
    }
    return lResult;
}

}  // namespace Integral
}  // namespace Math