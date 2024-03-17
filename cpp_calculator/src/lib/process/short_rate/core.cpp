/**
 * @file core.cpp
 * @brief This implements short rate path classes.
 * @author kakune
 * @date 1/29/2024
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "process/short_rate.hpp"

namespace Process
{
namespace ShortRate
{

void ModelAbstract::buildSpot()
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
}
void ModelAbstract::buildZCB()
{
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

void ModelAbstract::build()
{
    buildSpot();
    buildZCB();
}

double ModelAbstract::priceZCB( double inMaturityTime ) const
{
    return mInterpZCB( inMaturityTime );
}

double ModelAbstract::priceZCB( double inStartTime,
                                double inMaturityTime ) const
{
    return priceZCB( inMaturityTime ) / priceZCB( inStartTime );
}
double ModelAbstract::priceZCB( std::size_t inIndStartTime,
                                std::size_t inIndMaturityTime ) const
{
    if ( inIndStartTime >= msTerms->size() ||
         inIndMaturityTime >= msTerms->size() )
    {
        std::cerr << "Error: "
                     "Process::ShortRate::ModelAbstract:priceZCB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return priceZCB( msTerms->at( inIndStartTime ),
                     msTerms->at( inIndMaturityTime ) );
}

double ModelAbstract::forwardRate( double inStartTime,
                                   double inTerminalTime ) const
{
    return ( std::log( priceZCB( msTerms->at( 0 ), inStartTime ) ) -
             std::log( priceZCB( msTerms->at( 0 ), inTerminalTime ) ) ) /
           ( inTerminalTime - inStartTime );
}
double ModelAbstract::forwardRate( std::size_t inIndStartTime,
                                   std::size_t inIndTerminalTime ) const
{
    if ( inIndStartTime >= msTerms->size() ||
         inIndTerminalTime >= msTerms->size() )
    {
        std::cerr << "Error: "
                     "Process::ShortRate::ModelAbstract:forwardRate()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return forwardRate( msTerms->at( inIndStartTime ),
                        msTerms->at( inIndTerminalTime ) );
}

double ModelAbstract::instantaneousForwardRate( double inTime ) const
{
    return -mInterpZCB.deriv( inTime, 1 ) / mInterpZCB( inTime );
}
double ModelAbstract::instantaneousForwardRate( std::size_t inIndTime ) const
{
    if ( inIndTime >= msTerms->size() )
    {
        std::cerr
            << "Error: "
               "Process::ShortRate::ModelAbstract::instantaneousForwardRate()"
            << std::endl
            << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return instantaneousForwardRate( msTerms->at( inIndTime ) );
}

double ConstantRate::driftCoeff( std::size_t inIndPath,
                                 std::size_t inIndTerm ) const
{
    return 0.0;
}

void OneFactorAbstract::buildSpot()
{
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        muRandomPath->setIndexTime( iTerm );
        double lTmpDt     = msTerms->at( iTerm ) - msTerms->at( iTerm - 1 );
        double lHalfTmpDt = 0.5 * lTmpDt;
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mSpotRates.at( iPath ).at( iTerm ) =
                mSpotRates.at( iPath ).at( iTerm - 1 ) +
                driftCoeff( iPath, iTerm ) * lTmpDt +
                muRandomPath->generateRandomVal() * volCoeff( iPath, iTerm );
            mDFs.at( iPath ).at( iTerm ) =
                mDFs.at( iPath ).at( iTerm - 1 ) *
                std::exp( -( mSpotRates.at( iPath ).at( iTerm - 1 ) +
                             mSpotRates.at( iPath ).at( iTerm ) ) *
                          lHalfTmpDt );
        }
    }
}

}  // namespace ShortRate
}  // namespace Process