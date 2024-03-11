/**
 * @file one-factor_Gauss.cpp
 * @brief This implements one-factor Gaussian short-rate models. (see Chap.10 of
 * Andersen & Piterbarg)
 * @author kakune
 * @date 3/8/2024
 */

#include "process/short_rate/one-factor_Gauss.hpp"

#include <iostream>
#include <memory>

namespace Process
{
namespace ShortRate
{

double HoLee::driftCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    if ( msMarketData == nullptr ) { return 0.0; }
    return msMarketData->mInterpInstantaneousForwardRate.deriv(
               msTerms->at( inIndTerm - 1 ), 1 ) +
           mVol2 * ( msTerms->at( inIndTerm - 1 ) - msTerms->at( 0 ) );
}
double HoLee::volCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    return mVol;
}

double HoLee::analyticalPriceZCB( double inStartTime,
                                  double inTerminalTime ) const
{
    if ( msMarketData != nullptr )
    {
        return msMarketData->mInterpZCB( inTerminalTime ) /
               msMarketData->mInterpZCB( inStartTime );
    }
    double lTimeDif = inTerminalTime - inStartTime;
    return std::exp( -mInitSpotRate * lTimeDif +
                     mVol2 * lTimeDif * lTimeDif * lTimeDif / 6.0 );
}
double HoLee::analyticalPriceZCB( std::size_t inIndStartTime,
                                  std::size_t inIndMaturityTime ) const
{
    if ( inIndStartTime >= msTerms->size() ||
         inIndMaturityTime >= msTerms->size() )
    {
        std::cerr << "Error: "
                     "Process::ShortRate::HoLee:analyticalPriceZCB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
    }
    return analyticalPriceZCB( msTerms->at( inIndStartTime ),
                               msTerms->at( inIndMaturityTime ) );
}

double Vasicek::driftCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    --inIndTerm;
    if ( msMarketData == nullptr )
    {
        return mKappa * ( mMean - mSpotRates.at( inIndPath ).at( inIndTerm ) );
    }
    double lTime = msTerms->at( inIndTerm );
    double lMean =
        mMeanFactor1 *
            msMarketData->mInterpInstantaneousForwardRate.deriv( lTime, 1 ) +
        msMarketData->mInterpInstantaneousForwardRate( lTime ) +
        mMeanFactor2 *
            ( 1.0 - std::exp( -2.0 * mKappa * ( lTime - msTerms->at( 0 ) ) ) );
    return mKappa * ( lMean - mSpotRates.at( inIndPath ).at( inIndTerm ) );
}
double Vasicek::volCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    return mVol;
}

double Vasicek::analyticalPriceZCB( double inTerminalTime ) const
{
    if ( msMarketData != nullptr )
    {
        return msMarketData->mInterpZCB( inTerminalTime );
    }
    double lTimeDif = inTerminalTime - msTerms->at( 0 );
    double lB       = ( 1.0 - std::exp( -mKappa * lTimeDif ) ) / mKappa;
    double lA       = ( mMean - mVol * mVol / ( 2.0 * mKappa * mKappa ) ) *
                    ( lB - lTimeDif ) -
                ( mVol * mVol * lB * lB ) / ( 4.0 * mKappa );
    return std::exp( lA - lB * mInitSpotRate );
}
double Vasicek::analyticalPriceZCB( double inStartTime,
                                    double inTerminalTime ) const
{
    return analyticalPriceZCB( inTerminalTime ) /
           analyticalPriceZCB( inStartTime );
}
double Vasicek::analyticalPriceZCB( std::size_t inIndStartTime,
                                    std::size_t inIndMaturityTime ) const
{
    if ( inIndStartTime >= msTerms->size() ||
         inIndMaturityTime >= msTerms->size() )
    {
        std::cerr << "Error: "
                     "Process::ShortRate::HoLee:analyticalPriceZCB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
    }
    return analyticalPriceZCB( msTerms->at( inIndStartTime ),
                               msTerms->at( inIndMaturityTime ) );
}

void GSR::buildSpot()
{
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        muRandomPath->setIndexTime( iTerm );
        double lTmpDt     = msTerms->at( iTerm ) - msTerms->at( iTerm - 1 );
        double lHalfTmpDt = 0.5 * lTmpDt;
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            if ( msMarketData != nullptr )
            {
                mFactors.at( iPath ).at( iTerm ) =
                    mFactors.at( iPath ).at( iTerm - 1 ) +
                    factorCoeff( iPath, iTerm );
            }
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

double GSR::driftCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    --inIndTerm;
    double lTime = msTerms->at( inIndTerm );
    if ( msMarketData == nullptr )
    {
        return mInterpKappa( lTime ) *
               ( mInterpMean( lTime ) -
                 mSpotRates.at( inIndPath ).at( inIndTerm ) );
    }
    double lKappa = mInterpKappa( lTime );
    double lMean =
        msMarketData->mInterpInstantaneousForwardRate.deriv( lTime, 1 ) /
            lKappa +
        msMarketData->mInterpInstantaneousForwardRate( lTime ) +
        mFactors.at( inIndPath ).at( inIndTerm ) / lKappa;

    return lKappa * ( lMean - mSpotRates.at( inIndPath ).at( inIndTerm ) );
}
double GSR::volCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    return mInterpVol( msTerms->at( inIndTerm - 1 ) );
}

double GSR::factorCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    double lTime        = msTerms->at( inIndTerm - 1 );
    double lDt          = msTerms->at( inIndTerm ) - lTime;
    double lCut         = 0.05 * lDt;
    double lTmpKappa    = 0.5 * ( mInterpKappa( lTime + lCut ) +
                               mInterpKappa( lTime + lDt - lCut ) );
    double lTmpInvKappa = 1.0 / lTmpKappa;
    double lVol         = mInterpVol( lTime + lCut );
    double lDVol        = mInterpVol.deriv( lTime + lCut, 1 );
    double lResult      = 0.5 * ( lVol - lTmpInvKappa * lDVol ) *
                     ( 1.0 - std::exp( -2.0 * lTmpKappa * lDt ) );
    lResult += lDt * lDVol;
    return lResult * lVol * lTmpInvKappa;
}

}  // namespace ShortRate
}  // namespace Process
