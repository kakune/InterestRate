/**
 * @file core.cpp
 * @brief This implements short rate path classes.
 * @author kakune
 * @date 1/29/2024
 */

#include "short_rate/one-factor/core.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace ShortRate
{
namespace OneFactor
{

Process::MarketData::SpotRates ModelAbstract::calcSpotRates() const
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
    return Process::MarketData::SpotRates( mTerms, lSpots );
}

double ConstantRate::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return 0.0;
}

Process::MarketData::SpotRates OneFactorAbstract::calcSpotRates() const
{
    std::vector<std::vector<double>> lSpots(
        mNPath, std::vector<double>( mTerms.size(), mInitSpotRate ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpDt = mTerms[iTerm] - mTerms[iTerm - 1];
        muRandomPath->setIndexTime( iTerm );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lSpots[iPath][iTerm] = lSpots[iPath][iTerm - 1] +
                                   driftCoeff( iPath, iTerm, lSpots ) * lTmpDt +
                                   muRandomPath->generateRandomVal() *
                                       volCoeff( iPath, iTerm, lSpots );
        }
    }
    return Process::MarketData::SpotRates( mTerms, lSpots );
}

}  // namespace OneFactor
}  // namespace ShortRate