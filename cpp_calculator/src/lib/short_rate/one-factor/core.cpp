/**
 * @file core.cpp
 * @brief This implements short rate path classes.
 * @author kakune
 * @date 1/29/2024
 */

#include "short_rate/one-factor/core.hpp"

#include <cmath>
#include <vector>

namespace ShortRate
{
namespace OneFactor
{
#ifndef USE_CUDA
ShortRate::SpotRates ModelAbstract::createSpotRates() const
{
    std::vector<std::vector<double>> lSpots(
        mNPath, std::vector<double>( mTerms.size(), mInitSpotRate ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpDt = mTerms[iTerm] - mTerms[iTerm - 1];
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lSpots[iPath][iTerm] = lSpots[iPath][iTerm - 1] +
                                   driftCoeff( iPath, iTerm, lSpots ) * lTmpDt;
        }
    }
    return ShortRate::SpotRates( mTerms, lSpots );
}

double ConstantRate::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return 0.0;
}

ShortRate::SpotRates OneFactorAbstract::createSpotRates() const
{
    std::vector<std::vector<double>> lSpots(
        mNPath, std::vector<double>( mTerms.size(), mInitSpotRate ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpDt = mTerms[iTerm] - mTerms[iTerm - 1];
        muStdBrown->initialize();
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lSpots[iPath][iTerm] = lSpots[iPath][iTerm - 1] +
                                   driftCoeff( iPath, iTerm, lSpots ) * lTmpDt +
                                   ( *muStdBrown )() *
                                       mTerms.sqrtDifTime( iTerm ) *
                                       volCoeff( iPath, iTerm, lSpots );
        }
    }
    return ShortRate::SpotRates( mTerms, lSpots );
}
#endif
}  // namespace OneFactor
}  // namespace ShortRate