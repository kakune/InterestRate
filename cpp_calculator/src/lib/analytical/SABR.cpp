/**
 * @file SABR.hpp
 * @brief This implements analytical calculators relating to SABR Model.
 * @author kakune
 * @date 4/30/2024
 */

#include "analytical/SABR.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "analytical/Black76.hpp"
#include "math/findroot_1d.hpp"
#include "math/integral_1d.hpp"
#include "math/interpolate_1d.hpp"
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

static constexpr double gSqrt8      = 2.82842712475;
static constexpr double gDoublePI   = 2.0 * M_PI;
static constexpr double gMaxKernelS = 20.0;
static double kernelG( double inTau, double inS )
{
    if ( inS > gMaxKernelS ) { return 0.0; }
    const double lFactor = gSqrt8 * exp( -0.125 * inTau ) /
                           ( inTau * std::sqrt( gDoublePI * inTau ) );
    const double lCoshS        = cosh( inS );
    const double lInvDoubleTau = 0.5 / inTau;
    auto lIntegrand            = [lCoshS, lInvDoubleTau]( double inU ) -> double
    {
        if ( inU > gMaxKernelS ) { return 0.0; }
        return inU * exp( -inU * inU * lInvDoubleTau ) *
               std::sqrt( std::max( { 0.0, cosh( inU ) - lCoshS } ) );
    };
    return lFactor * Math::Integral::UpperInfiniteInterval::DEFormulaForExp(
                         lIntegrand, inS, 0.001 );
}

static auto makeKernelG( double inTau, double inMin, double inEpsRel = 1e-8,
                         std::size_t inNDeg = 3 )
{
    double lDx = std::sqrt( std::abs( inTau ) ) / 200.0;
    std::vector<double> lXs;
    for ( double lX = inMin - 0.001; lX < gMaxKernelS; lX += lDx )
    {
        lXs.emplace_back( lX );
    }
    lXs.emplace_back( gMaxKernelS + 0.001 );
    std::vector<double> lYs( lXs.size(), 0.0 );

    double lMaxX = gMaxKernelS;
    for ( std::size_t iX = 0; iX < lXs.size() - 1; ++iX )
    {
        lYs[iX] = kernelG( inTau, lXs[iX] );
        if ( lYs[iX] < lYs[0] * inEpsRel )
        {
            lMaxX = lXs[iX + 1];
            break;
        }
    }
    Math::Interpolate1D::NewtonSpline lSpline(
        std::make_shared<std::vector<double>>( lXs ),
        std::make_shared<std::vector<double>>( lYs ), inNDeg );

    auto lResult = [lSpline, lMaxX]( double inS )
    {
        if ( inS >= lMaxX ) { return 0.0; }
        return lSpline( inS );
    };
    return lResult;
}

double OneTerm::approxBlackImpVolByAntonov() const
{
    if ( mStrike == mInitPrice )
    {
        OneTerm lTmpModelPlus  = *this;
        OneTerm lTmpModelMinus = *this;
        double lDStrike        = mStrike * 1e-2;
        lTmpModelMinus.setStrike( mStrike - lDStrike );
        lTmpModelPlus.setStrike( mStrike + lDStrike );
        return 0.5 * ( lTmpModelMinus.approxBlackImpVolByAntonov() +
                       lTmpModelPlus.approxBlackImpVolByAntonov() );
    }
    const double lOneMinusExponent    = 1.0 - mExponent;
    const double lInvOneMinusExponent = 1.0 / lOneMinusExponent;
    const double lPowInitPrice = std::pow( mInitPrice, lOneMinusExponent );
    const double lPowStrike    = std::pow( mStrike, lOneMinusExponent );
    const double lSquareTildeVolVol =
        ( mVolVol - 1.5 * mCorr *
                        ( mVolVol * mCorr +
                          mInitVol * lOneMinusExponent / lPowInitPrice ) );
    if ( lSquareTildeVolVol <= 0.0 )
    {
        throw std::invalid_argument(
            std::string(
                "Analytical::SABR::OneTerm::approxBlackImpVolByAntonov()\n" ) +
            std::string( "Correlation is too large." ) );
    }
    const double lTildeVolVol = std::sqrt( lSquareTildeVolVol );
    const double lQ           = lPowStrike * lInvOneMinusExponent;
    const double lQZero       = lPowInitPrice * lInvOneMinusExponent;
    const double lEta         = std::abs( 0.5 * lInvOneMinusExponent );

    const double lDeltaQ       = lQ - lQZero;
    const double lVolVolDeltaQ = mVolVol * lDeltaQ;
    const double lMinInitVol =
        std::sqrt( lVolVolDeltaQ * ( lVolVolDeltaQ + 2.0 * mCorr * mInitVol ) +
                   mInitVol * mInitVol );
    const double lPhi =
        std::pow( ( lMinInitVol + mCorr * mInitVol + mVolVol * lDeltaQ ) /
                      ( ( 1.0 + mCorr ) * mInitVol ),
                  lTildeVolVol / mVolVol );
    const double lTildeAlphaZero =
        2.0 * lPhi * lDeltaQ * lTildeVolVol / ( lPhi * lPhi - 1.0 );

    const double lSqrtCorr = std::sqrt( 1.0 - mCorr * mCorr );
    const double lUZero =
        ( lDeltaQ * mVolVol * mCorr + mInitVol - lMinInitVol ) /
        ( lDeltaQ * mVolVol * lSqrtCorr );
    const double lL =
        lMinInitVol * lOneMinusExponent / ( lPowStrike * mVolVol * lSqrtCorr );
    const double lInvSqrtL = 1.0 / std::sqrt( std::abs( 1.0 - lL * lL ) );
    const double lI =
        ( lL < 1 )
            ? 2.0 * lInvSqrtL *
                  ( atan( ( lUZero + lL ) * lInvSqrtL ) -
                    atan( lL * lInvSqrtL ) )
            : lInvSqrtL * log( ( ( 1.0 + lUZero * lL ) * lInvSqrtL + lUZero ) /
                               ( ( 1.0 + lUZero * lL ) * lInvSqrtL - lUZero ) );
    const double lPhiZero =
        acos( -( lDeltaQ * mVolVol + mInitVol * mCorr ) / lMinInitVol );
    const double lMinB =
        -0.5 * ( ( mExponent * mCorr ) / ( lOneMinusExponent * lSqrtCorr ) ) *
        ( M_PI - lPhiZero - acos( mCorr ) - lI );

    const double lTildeAlphaOne =
        ( 0.5 * lTildeAlphaZero * lTildeVolVol * lTildeVolVol /
          ( log( lPhi ) * ( lPhi * lPhi - 1.0 ) / ( lPhi * lPhi + 1.0 ) ) ) *
        ( log( mInitVol * lMinInitVol ) -
          log( lTildeAlphaZero *
               std::sqrt( lDeltaQ * lDeltaQ * lTildeVolVol * lTildeVolVol +
                          lTildeAlphaZero * lTildeAlphaZero ) ) -
          2.0 * lMinB );

    const double lTildeAlpha    = lTildeAlphaZero + mTime * lTildeAlphaOne;
    const double lInvTildeAlpha = 1.0 / lTildeAlpha;
    const double lSMinus =
        asinh( lTildeVolVol * std::abs( lDeltaQ ) * lInvTildeAlpha );
    const double lSPlus =
        asinh( lTildeVolVol * std::abs( lQ + lQZero ) * lInvTildeAlpha );
    const double lSquareSinhSMinus = sinh( lSMinus ) * sinh( lSMinus );
    const double lSquareSinhSPlus  = sinh( lSPlus ) * sinh( lSPlus );

    const auto lKernelG =
        makeKernelG( mTime * lTildeVolVol * lTildeVolVol, lSMinus );
    // std::cout << "Kernel created." << std::endl;
    auto lIntegrand1 = [&]( double inS ) -> double
    {
        const double lG = lKernelG( inS );
        if ( lG == 0.0 ) { return 0.0; }
        const double lSinh       = sinh( inS );
        const double lSquareSinh = lSinh * lSinh;
        const double lKappa =
            2.0 * atan( std::sqrt( std::max(
                      { 0.0, ( lSquareSinh - lSquareSinhSMinus ) /
                                 ( lSquareSinhSPlus - lSquareSinh ) } ) ) );
        return lG * sin( lEta * lKappa ) / lSinh;
    };

    const double lSinEtaPi = sin( lEta * M_PI );
    auto lIntegrand2       = [&]( double inS ) -> double
    {
        const double lG = lKernelG( inS );
        if ( lG == 0.0 ) { return 0.0; }
        const double lSinh       = sinh( inS );
        const double lSquareSinh = lSinh * lSinh;
        const double lPsi =
            2.0 * atanh( std::sqrt( ( lSquareSinh - lSquareSinhSPlus ) /
                                    ( lSquareSinh - lSquareSinhSMinus ) ) );
        return lSinEtaPi * std::exp( -lEta * lPsi ) * lG / lSinh;
    };

    const double lPriceBlack =
        std::max( { 0.0, mInitPrice - mStrike } ) +
        M_2_PI * std::sqrt( mInitPrice * mStrike ) *
            ( Math::Integral::FiniteInterval::DEFormula( lIntegrand1, lSMinus,
                                                         lSPlus, 0.001 ) +
              Math::Integral::UpperInfiniteInterval::DEFormula(
                  lIntegrand2, lSPlus, 0.001 ) );
    auto lRootGrand = [this, lPriceBlack]( double inVol )
    {
        return Black76::funcBlackPositive( mStrike, mInitPrice,
                                           inVol * std::sqrt( mTime ) ) -
               lPriceBlack;
    };
    return Math::FindRoot1D::Brent( lRootGrand, 1e-10, 100 );
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