/**
 * @file one-factor_Gauss.cpp
 * @brief This implements one-factor Gaussian short-rate models. (see Chap.10 of
 * Andersen & Piterbarg)
 * @author kakune
 * @date 3/8/2024
 */

#include <iostream>
#include <memory>

#include "process/short_rate_MC.hpp"

namespace Process
{
namespace ShortRateMC
{

double HoLee::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return 0.0;
    // return msMarketData->mInterpInstantaneousForwardRate.deriv(
    //            mTerms[inIndTerm - 1], 1 ) +
    //        mVol2 * ( mTerms[inIndTerm - 1] - mTerms.front() );
}
double HoLee::volCoeff( std::size_t inIndPath, std::size_t inIndTerm,
                        const std::vector<std::vector<double>>& inSpots ) const
{
    return mVol;
}

double HoLee::analyticalPriceZCB( double inStartTime,
                                  double inTerminalTime ) const
{
    double lTimeDif = inTerminalTime - inStartTime;
    return std::exp( -mInitSpotRate * lTimeDif +
                     mVol * mVol * lTimeDif * lTimeDif * lTimeDif / 6.0 );
}
double HoLee::analyticalPriceZCB( std::size_t inIndStartTime,
                                  std::size_t inIndMaturityTime ) const
{
    if ( inIndStartTime >= mTerms.size() || inIndMaturityTime >= mTerms.size() )
    {
        std::cerr << "Error: "
                     "Process::ShortRate::HoLee:analyticalPriceZCB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return analyticalPriceZCB( mTerms[inIndStartTime],
                               mTerms[inIndMaturityTime] );
}

double HoLeeWithMarket::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return mMarketZCB.derivInstantaneousForwardRate( mTerms[inIndTerm - 1] ) +
           mVol2 * ( mTerms[inIndTerm - 1] - mTerms.front() );
}
double HoLeeWithMarket::volCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return mVol;
}

double Vasicek::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return mKappa * ( mMean - inSpots[inIndPath][inIndTerm - 1] );
}
double Vasicek::volCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return mVol;
}

double Vasicek::analyticalPriceZCB( double inTerminalTime ) const
{
    double lTimeDif = inTerminalTime - mTerms.front();
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
    if ( inIndStartTime >= mTerms.size() || inIndMaturityTime >= mTerms.size() )
    {
        std::cerr << "Error: "
                     "Process::ShortRate::HoLee:analyticalPriceZCB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return analyticalPriceZCB( mTerms[inIndStartTime],
                               mTerms[inIndMaturityTime] );
}

double VasicekWithMarket::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    --inIndTerm;
    double lTime = mTerms[inIndTerm];
    double lMean =
        mMeanFactor1 * mMarketZCB.derivInstantaneousForwardRate( lTime ) +
        mMarketZCB.instantaneousForwardRate( lTime ) +
        mMeanFactor2 *
            ( 1.0 - std::exp( -2.0 * mKappa * ( lTime - mTerms.front() ) ) );
    return mKappa * ( lMean - inSpots[inIndPath][inIndTerm] );
}
double VasicekWithMarket::volCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return mVol;
}

double GSR::driftCoeff( std::size_t inIndPath, std::size_t inIndTerm,
                        const std::vector<std::vector<double>>& inSpots ) const
{
    --inIndTerm;
    double lTime = mTerms[inIndTerm];
    return mInterpKappa( lTime ) *
           ( mInterpMean( lTime ) - inSpots[inIndPath][inIndTerm] );
}
double GSR::volCoeff( std::size_t inIndPath, std::size_t inIndTerm,
                      const std::vector<std::vector<double>>& inSpots ) const
{
    return mInterpVol( mTerms[inIndTerm - 1] );
}

MarketData::SpotRates GSRWithMarket::calcSpotRates() const
{
    std::vector<std::vector<double>> lSpots(
        mNPath, std::vector<double>( mTerms.size(), mInitSpotRate ) );
    std::vector<std::vector<double>> lFactors(
        mNPath, std::vector<double>( mTerms.size(), 0.0 ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        double lTmpDt = mTerms[iTerm] - mTerms[iTerm - 1];
        muRandomPath->setIndexTime( iTerm );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lFactors[iPath][iTerm] =
                lFactors[iPath][iTerm - 1] +
                factorCoeff( iPath, iTerm, lSpots, lFactors );
            lSpots[iPath][iTerm] =
                lSpots[iPath][iTerm - 1] +
                driftCoeff( iPath, iTerm, lSpots, lFactors ) * lTmpDt +
                muRandomPath->generateRandomVal() *
                    volCoeff( iPath, iTerm, lSpots, lFactors );
        }
    }
    return MarketData::SpotRates( mTerms, lSpots );
}

double GSRWithMarket::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots,
    const std::vector<std::vector<double>>& inFactors ) const
{
    --inIndTerm;
    double lTime  = mTerms[inIndTerm];
    double lKappa = mInterpKappa( lTime );
    double lMean  = mMarketZCB.derivInstantaneousForwardRate( lTime ) / lKappa +
                   mMarketZCB.instantaneousForwardRate( lTime ) +
                   inFactors[inIndPath][inIndTerm] / lKappa;

    return lKappa * ( lMean - inSpots[inIndPath][inIndTerm] );
}
double GSRWithMarket::volCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots,
    const std::vector<std::vector<double>>& inFactors ) const
{
    return mInterpVol( mTerms[inIndTerm - 1] );
}

double GSRWithMarket::factorCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots,
    const std::vector<std::vector<double>>& inFactors ) const
{
    double lTime        = mTerms[inIndTerm - 1];
    double lDt          = mTerms[inIndTerm] - lTime;
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

}  // namespace ShortRateMC
}  // namespace Process
