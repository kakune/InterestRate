/**
 * @file core.hpp
 * @brief This implements classes to calculate multi-factor short rate paths.
 * @author kakune
 * @date 3/23/2024
 */

#include "short_rate/multi-factor/core.hpp"

#include <cmath>
#include <vector>

namespace ShortRate
{
namespace MultiFactor
{

ShortRate::SpotRates ModelAbstract::createSpotRates() const
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
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lRates[iPath][iTerm] =
                transfStateToRate( lSpots[iPath][iTerm], iTerm );
        }
    }
    return ShortRate::SpotRates( mTerms, lRates );
}
double ModelAbstract::transfStateToRate( const Math::Vec& inState,
                                         std::size_t inIndTime ) const
{
    return inState[0];
}

Math::Vec ConstantRate::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots ) const
{
    return Math::makeVec( mDim, 0.0 );
}

ShortRate::SpotRates MultiFactorAbstract::createSpotRates() const
{
    std::vector<std::vector<Math::Vec>> lSpots(
        mNPath, std::vector<Math::Vec>( mTerms.size(), mInitState ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpDt = mTerms[iTerm] - mTerms[iTerm - 1];
        muStdBrown->initialize();
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lSpots[iPath][iTerm] =
                lSpots[iPath][iTerm - 1] +
                driftCoeff( iPath, iTerm, lSpots ) * lTmpDt +
                volTerm( iPath, iTerm, lSpots,
                         ( *muStdBrown )() * mTerms.sqrtDifTime( iTerm ) );
        }
    }
    std::vector<std::vector<double>> lRates(
        mNPath, std::vector<double>( mTerms.size() ) );
    for ( std::size_t iTerm = 0; iTerm < mTerms.size(); ++iTerm )
    {
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lRates[iPath][iTerm] =
                transfStateToRate( lSpots[iPath][iTerm], iTerm );
        }
    }
    return ShortRate::SpotRates( mTerms, lRates );
}

}  // namespace MultiFactor
}  // namespace ShortRate