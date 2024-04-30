/**
 * @file SABR.hpp
 * @brief This implements analytical calculators relating to SABR Model.
 * @author kakune
 * @date 4/30/2024
 */

#include "analytical/SABR.hpp"

#include <iostream>

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

}  // namespace Analytical::SABR