/**
 * @file SABR.hpp
 * @brief This implements analytical calculators relating to SABR Model.
 * @author kakune
 * @date 4/30/2024
 */

#include "analytical/SABR.hpp"

#include <iostream>
#include <stdexcept>

#include "math/matrix.hpp"
#include "math/optimize.hpp"

namespace Analytical::SABR
{

OneTerm::OneTerm( double inTime ) : mTime( inTime ) {}
OneTerm& OneTerm::setInitPrice( double inInitPrice )
{
    mInitPrice = inInitPrice;
    return *this;
}
OneTerm& OneTerm::setInitVol( double inInitVol )
{
    mInitVol = inInitVol;
    return *this;
}
OneTerm& OneTerm::setExponent( double inExponent )
{
    mExponent = inExponent;
    return *this;
}
OneTerm& OneTerm::setCorr( double inCorr )
{
    mCorr = inCorr;
    return *this;
}
OneTerm& OneTerm::setVolVol( double inVolVol )
{
    mVolVol = inVolVol;
    return *this;
}
OneTerm& OneTerm::setStrike( double inStrike )
{
    mStrike = inStrike;
    return *this;
}
const OneTerm& OneTerm::printAllParam() const
{
    std::cout << "InitVol  : " << mInitVol << std::endl;
    std::cout << "Exponent : " << mExponent << std::endl;
    std::cout << "Corr     : " << mCorr << std::endl;
    std::cout << "VolVol   : " << mVolVol << std::endl;
    return *this;
}
double OneTerm::approxBlackImpVolByHaganATM() const
{
    const double lOneMinusBeta = 1.0 - mExponent;
    const double lInvInitPricePow =
        std::pow( mInitPrice, -lOneMinusBeta );  // 1 / F^(1-beta)
    const double lTmp     = lOneMinusBeta * mInitVol * lInvInitPricePow;
    const double lFactor1 = gInv24 * lTmp * lTmp;
    const double lFactor2 =
        0.25 * mInitVol * mExponent * mVolVol * mCorr * lInvInitPricePow;
    const double lFactor3 =
        gInv24 * ( 2.0 - 3.0 * mCorr * mCorr ) * mVolVol * mVolVol;
    return mInitVol * lInvInitPricePow *
           ( 1.0 + ( lFactor1 + lFactor2 + lFactor3 ) * mTime );
}
double OneTerm::approxBlackImpVolByHagan() const
{
    if ( mStrike == mInitPrice ) { return approxBlackImpVolByHaganATM(); }
    const double lOneMinusBeta          = 1.0 - mExponent;
    const double lInvInitPriceStrikePow = std::pow(
        mStrike * mInitPrice, -0.5 * lOneMinusBeta );  // 1/(FK)^((1-beta)/2)
    const double lLogPriceStrike = std::log( mInitPrice / mStrike );
    const double lZ =
        ( mVolVol * lLogPriceStrike ) / ( mInitVol * lInvInitPriceStrikePow );
    const double lXi = std::log(
        ( std::sqrt( 1.0 + ( -2.0 * mCorr + lZ ) * lZ ) + lZ - mCorr ) /
        ( 1.0 - mCorr ) );
    if ( lXi == 0.0 ) { return approxBlackImpVolByHaganATM(); }
    const double lTmp     = lOneMinusBeta * mInitVol * lInvInitPriceStrikePow;
    const double lFactor1 = gInv24 * lTmp * lTmp;
    const double lFactor2 =
        0.25 * mInitVol * mExponent * mVolVol * mCorr * lInvInitPriceStrikePow;
    const double lFactor3 =
        gInv24 * ( 2.0 - 3.0 * mCorr * mCorr ) * mVolVol * mVolVol;
    const double lDenomFactor = ( lOneMinusBeta * lLogPriceStrike ) *
                                ( lOneMinusBeta * lLogPriceStrike );
    return mInitVol * lInvInitPriceStrikePow * lZ *
           ( 1.0 + ( lFactor1 + lFactor2 + lFactor3 ) * mTime ) /
           ( lXi *
             ( 1.0 + ( gInv24 + gInv1920 * lDenomFactor ) * lDenomFactor ) );
}

double OneTerm::approxNormalImpVolByHagan() const
{
    const double lOneMinusBeta = 1.0 - mExponent;
    const double lFactor3 =
        gInv24 * ( 2.0 - 3.0 * mCorr * mCorr ) * mVolVol * mVolVol;
    if ( mStrike == mInitPrice )
    {
        const double lInvInitPricePow =
            std::pow( mInitPrice, -lOneMinusBeta );  // 1 / F^(1-beta)
        const double lTmp = mInitVol * lInvInitPricePow;
        const double lFactor1 =
            -gInv24 * mExponent * ( 2.0 - mExponent ) * lTmp * lTmp;
        const double lFactor2 =
            0.25 * mInitVol * mExponent * mVolVol * mCorr * lInvInitPricePow;
        return mInitVol * ( mInitPrice / lInvInitPricePow ) *
               ( 1.0 + ( lFactor1 + lFactor2 + lFactor3 ) * mTime );
    }
    const double lPriceStrike = mStrike * mInitPrice;
    const double lInvInitPriceStrikePow =
        std::pow( lPriceStrike, -0.5 * lOneMinusBeta );  // 1/(FK)^((1-beta)/2)
    const double lLogPriceStrike = std::log( mInitPrice / mStrike );
    const double lZ =
        ( mVolVol * lLogPriceStrike ) / ( mInitVol * lInvInitPriceStrikePow );
    const double lXi = std::log(
        ( std::sqrt( 1.0 + ( -2.0 * mCorr + lZ ) * lZ ) + lZ - mCorr ) /
        ( 1.0 - mCorr ) );
    const double lTmp = mInitVol * lInvInitPriceStrikePow;
    const double lFactor1 =
        -gInv24 * mExponent * ( 2.0 - mExponent ) * lTmp * lTmp;
    const double lFactor2 =
        0.25 * mInitVol * mExponent * mVolVol * mCorr * lInvInitPriceStrikePow;
    const double lSquareLogPriceStrike = lLogPriceStrike * lLogPriceStrike;
    const double lDenomFactor =
        lOneMinusBeta * lOneMinusBeta * lSquareLogPriceStrike;

    return mInitVol * std::sqrt( lPriceStrike ) * lZ *
           ( 1.0 + ( lFactor1 + lFactor2 + lFactor3 ) * mTime ) *
           ( 1.0 + ( gInv24 + gInv1920 * lSquareLogPriceStrike ) *
                       lSquareLogPriceStrike ) /
           ( lXi * lInvInitPriceStrikePow *
             ( 1.0 + ( gInv24 + gInv1920 * lDenomFactor ) * lDenomFactor ) );
}

static std::vector<OneTerm> prepareModels(
    double inTime, const std::vector<double>& inInitPrices,
    const std::vector<double>& inStrikes ) noexcept
{
    std::size_t lNPoints = inInitPrices.size();
    std::vector<OneTerm> lModels( lNPoints, inTime );
    for ( std::size_t iModel = 0; iModel < lNPoints; ++iModel )
    {
        lModels[iModel].setInitPrice( inInitPrices[iModel] );
        lModels[iModel].setStrike( inStrikes[iModel] );
    }
    return lModels;
}

static double calcInitialInitVol( const std::vector<double>& inInitPrices,
                                  const std::vector<double>& inStrikes,
                                  const std::vector<double>& inImpVols,
                                  double inExponent )
{
    double lMinDif   = 1e18;
    std::size_t lInd = 0;
    for ( std::size_t i = 0; i < inInitPrices.size(); ++i )
    {
        double lTmpDif = std::abs( inInitPrices[i] - inStrikes[i] );
        if ( lTmpDif < lMinDif )
        {
            lMinDif = lTmpDif;
            lInd    = i;
        }
    }
    return inImpVols[lInd] * std::pow( inInitPrices[lInd], 1.0 - inExponent );
}

OneTerm calibrateAllParam( double inTime,
                           const std::vector<double>& inInitPrices,
                           const std::vector<double>& inStrikes,
                           const std::vector<double>& inImpVols,
                           std::vector<double> inInitialParam )
{
    std::size_t lNPoints = inInitPrices.size();
    if ( lNPoints != inStrikes.size() || lNPoints != inImpVols.size() )
    {
        throw std::invalid_argument(
            std::string( "Analytical::SABR::calibrateAllParam()\n" ) +
            std::string( "size of inputs do not match." ) );
    }
    std::vector<OneTerm> lModels =
        prepareModels( inTime, inInitPrices, inStrikes );
    Math::Vec lObjectiveImpVols = Math::makeVec( inImpVols );

    auto lFunc = [&lNPoints, &lModels,
                  &lObjectiveImpVols]( const Math::Vec& lParam ) -> double
    {
        double lResult = 0.0;
        for ( std::size_t iModel = 0; iModel < lNPoints; ++iModel )
        {
            lModels[iModel].setInitVol( lParam[0] );
            lModels[iModel].setExponent( lParam[1] );
            lModels[iModel].setCorr( lParam[2] );
            lModels[iModel].setVolVol( lParam[3] );
            double lDif = lModels[iModel].approxBlackImpVolByHagan() -
                          lObjectiveImpVols[iModel];
            lResult += lDif * lDif;
        }
        return lResult;
    };

    if ( inInitialParam[0] <= 0.0 )
    {
        inInitialParam[0] = calcInitialInitVol( inInitPrices, inStrikes,
                                                inImpVols, inInitialParam[1] );
    }
    Math::Vec lResultParam =
        Math::Optimize::adam( lFunc, Math::makeVec( inInitialParam ), 2e-6,
                              0.85, 0.999, 1e-12, 1e7, 1e-20 );

    OneTerm lResultModel( inTime );
    lResultModel.setInitVol( lResultParam[0] );
    lResultModel.setExponent( lResultParam[1] );
    lResultModel.setCorr( lResultParam[2] );
    lResultModel.setVolVol( lResultParam[3] );

    return lResultModel;
}

OneTerm calibrateAllParam( double inTime, double inInitPrice,
                           const std::vector<double>& inStrikes,
                           const std::vector<double>& inImpVols,
                           std::vector<double> inInitialParam )
{
    return calibrateAllParam(
        inTime, std::vector<double>( inStrikes.size(), inInitPrice ), inStrikes,
        inImpVols, inInitialParam );
}

OneTerm calibrateParamWithFixedExponent(
    double inTime, const std::vector<double>& inInitPrices,
    const std::vector<double>& inStrikes, const std::vector<double>& inImpVols,
    double inExponent, std::vector<double> inInitialParam )
{
    std::size_t lNPoints = inInitPrices.size();
    if ( lNPoints != inStrikes.size() || lNPoints != inImpVols.size() )
    {
        throw std::invalid_argument(
            std::string( "Analytical::SABR::calibrateAllParam()\n" ) +
            std::string( "size of inputs do not match." ) );
    }
    std::vector<OneTerm> lModels =
        prepareModels( inTime, inInitPrices, inStrikes );
    for ( OneTerm& lModel : lModels ) { lModel.setExponent( inExponent ); }

    Math::Vec lObjectiveImpVols = Math::makeVec( inImpVols );

    auto lFunc = [&lNPoints, &lModels,
                  &lObjectiveImpVols]( const Math::Vec& lParam ) -> double
    {
        double lResult = 0.0;
        for ( std::size_t iModel = 0; iModel < lNPoints; ++iModel )
        {
            lModels[iModel].setInitVol( lParam[0] );
            lModels[iModel].setCorr( lParam[1] );
            lModels[iModel].setVolVol( lParam[2] );
            double lDif = lModels[iModel].approxBlackImpVolByHagan() -
                          lObjectiveImpVols[iModel];
            lResult += lDif * lDif;
        }
        return lResult;
    };

    if ( inInitialParam[0] <= 0.0 )
    {
        inInitialParam[0] = calcInitialInitVol( inInitPrices, inStrikes,
                                                inImpVols, inExponent );
    }
    Math::Vec lResultParam =
        Math::Optimize::adam( lFunc, Math::makeVec( inInitialParam ), 2e-6,
                              0.85, 0.999, 1e-12, 1e7, 1e-20 );

    OneTerm lResultModel( inTime );
    lResultModel.setInitVol( lResultParam[0] );
    lResultModel.setExponent( inExponent );
    lResultModel.setCorr( lResultParam[1] );
    lResultModel.setVolVol( lResultParam[2] );

    return lResultModel;
}

OneTerm calibrateParamWithFixedExponent( double inTime, double inInitPrice,
                                         const std::vector<double>& inStrikes,
                                         const std::vector<double>& inImpVols,
                                         double inExponent,
                                         std::vector<double> inInitialParam )
{
    return calibrateParamWithFixedExponent(
        inTime, std::vector<double>( inStrikes.size(), inInitPrice ), inStrikes,
        inImpVols, inExponent, inInitialParam );
}

}  // namespace Analytical::SABR