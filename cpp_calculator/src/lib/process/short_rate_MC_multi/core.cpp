/**
 * @file core.hpp
 * @brief This implements classes to calculate multi-factor short rate paths.
 * @author kakune
 * @date 3/23/2024
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "process/short_rate_MC.hpp"

namespace Process
{
namespace ShortRateMCMulti
{

MarketData::SpotRates ModelAbstract::calcSpotRates() const
{
    std::vector<std::vector<Math::Vec>> lSpots(
        mNPath, std::vector<Math::Vec>( mTerms.size(), mInitState ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpDt = mTerms[iTerm] - mTerms[iTerm - 1];
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lSpots[iPath][iTerm] = lSpots[iPath][iTerm - 1] +
                                   driftCoeff( iPath, iTerm, lSpots ) * lTmpDt;
        }
    }
    std::vector<std::vector<double>> lRates(
        mNPath, std::vector<double>( mTerms.size() ) );
    for ( std::size_t iTerm = 0; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpTime = mTerms[iTerm];
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lRates[iPath][iTerm] =
                transfStateToRate( lSpots[iPath][iTerm], lTmpTime );
        }
    }
    return MarketData::SpotRates( mTerms, lRates );
}
double ModelAbstract::transfStateToRate( const Math::Vec& inState,
                                         double inTime ) const
{
    return inState( 0 );
}

Math::Vec ConstantRate::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots ) const
{
    return Math::Vec( mDim, 0.0 );
}

MarketData::SpotRates MultiFactorAbstract::calcSpotRates() const
{
    std::vector<std::vector<Math::Vec>> lSpots(
        mNPath, std::vector<Math::Vec>( mTerms.size(), mInitState ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpDt = mTerms[iTerm] - mTerms[iTerm - 1];
        muRandomPath->setIndexTime( iTerm );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lSpots[iPath][iTerm] = lSpots[iPath][iTerm - 1] +
                                   driftCoeff( iPath, iTerm, lSpots ) * lTmpDt +
                                   volTerm( iPath, iTerm, lSpots,
                                            muRandomPath->generateRandomVal() );
        }
    }
    std::vector<std::vector<double>> lRates(
        mNPath, std::vector<double>( mTerms.size() ) );
    for ( std::size_t iTerm = 0; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpTime = mTerms[iTerm];
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lRates[iPath][iTerm] =
                transfStateToRate( lSpots[iPath][iTerm], lTmpTime );
        }
    }
    return MarketData::SpotRates( mTerms, lRates );
}

}  // namespace ShortRateMCMulti
}  // namespace Process