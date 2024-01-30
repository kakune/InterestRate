/**
 * @file short_rate.cpp
 * @brief This implements short rate path classes.
 * @author kakune
 * @date 1/29/2024
 */

#include "process/short_rate.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace Process
{
namespace ShortRate
{

void ModelAbstract::calcZCB()
{
    calcEachRates();
    std::vector< double > lIntegralRates( mNPath, 0 );
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        double lHalfTime =
            0.5 * ( msTerms->at( iTerm ) - msTerms->at( iTerm - 1 ) );
        mZCB[iTerm] = 0.0;
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lIntegralRates[iTerm] =
                lIntegralRates[iTerm - 1] +
                lHalfTime *
                    ( mSpotRates[iPath][iTerm] + mSpotRates[iPath][iTerm - 1] );
            mZCB[iTerm] += std::exp( -lIntegralRates[iTerm] );
        }
        mZCB[iTerm] /= mNPath;
    }
    mInterpZCB.build( msTerms,
                      std::make_shared< std::vector< double > >( mZCB ) );
}

double ModelAbstract::priceZCB( double inStartTime, double inMaturityTime )
{
    return mInterpZCB( inMaturityTime ) / mInterpZCB( inStartTime );
}

double ModelAbstract::forwardRate( double inStartTime, double inTerminalTime )
{
    auto lItr =
        std::upper_bound( msTerms->begin(), msTerms->end(), inTerminalTime );
    double lDTime = ( *lItr - inTerminalTime ) * 0.1;
    return ( std::log( priceZCB( inStartTime, inTerminalTime ) ) -
             std::log( priceZCB( inStartTime, inTerminalTime + lDTime ) ) ) /
           lDTime;
}

double ConstantRate::priceZCB( double inStartTime, double inMaturityTime )
{
    return std::exp( -( inMaturityTime - inStartTime ) * mRate );
}
double ConstantRate::forwardRate( double inStartTime, double inTerminalTime )
{
    return mRate;
}

}  // namespace ShortRate
}  // namespace Process