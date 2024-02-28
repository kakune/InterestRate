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

void ModelAbstract::build()
{
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        double lTmpDt     = msTerms->at( iTerm ) - msTerms->at( iTerm - 1 );
        double lHalfTmpDt = 0.5 * lTmpDt;
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mSpotRates.at( iPath ).at( iTerm ) =
                mSpotRates.at( iPath ).at( iTerm - 1 ) +
                driftCoeff( iPath, iTerm ) * lTmpDt;
            mDFs.at( iPath ).at( iTerm ) =
                mDFs.at( iPath ).at( iTerm - 1 ) *
                std::exp( -( mSpotRates.at( iPath ).at( iTerm - 1 ) +
                             mSpotRates.at( iPath ).at( iTerm ) ) *
                          lHalfTmpDt );
        }
    }

    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        mExpectedSpotRates.at( iTerm ) = 0.0;
        mZCBs.at( iTerm )              = 0.0;
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mExpectedSpotRates.at( iTerm ) +=
                mSpotRates.at( iPath ).at( iTerm );
            mZCBs.at( iTerm ) += mDFs.at( iPath ).at( iTerm );
        }
        mExpectedSpotRates.at( iTerm ) /= double( mNPath );
        mZCBs.at( iTerm ) /= double( mNPath );
    }
    mInterpZCB.build( msTerms,
                      std::make_shared<std::vector<double> >( mZCBs ) );
}

double ModelAbstract::priceZCB( double inStartTime, double inMaturityTime )
{
    return mInterpZCB( inMaturityTime ) / mInterpZCB( inStartTime );
}

double ModelAbstract::forwardRate( double inStartTime, double inTerminalTime )
{
    return ( std::log( priceZCB( msTerms->at( 0 ), inTerminalTime ) ) -
             std::log( priceZCB( msTerms->at( 0 ), inStartTime ) ) ) /
           ( inTerminalTime - inStartTime );
}

double ModelAbstract::instantaneousForwardRate( double inTime )
{
    return mInterpZCB.deriv( inTime, 1 ) / mInterpZCB( inTime );
}

double ConstantRate::driftCoeff( std::size_t inIndPath, std::size_t inIndTerm )
{
    return 0.0;
}

}  // namespace ShortRate
}  // namespace Process