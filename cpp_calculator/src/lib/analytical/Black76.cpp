/**
 * @file Black76.cpp
 * @brief This implements analytical calculators relating to Black-76 Model.
 * Model.
 * @author kakune
 * @date 4/13/2024
 */

#include "analytical/Black76.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "math/special_functions.hpp"

namespace Analytical
{
namespace Black76
{

double funcBlackPositive( double inStrike, double inPrice, double inVol )
{
    if ( inVol == 0.0 ) { return std::max( { inPrice - inStrike, 0.0 } ); }
    double lDPlus  = std::log( inPrice / inStrike ) / inVol + 0.5 * inVol;
    double lDMinus = lDPlus - inVol;
    return ( inPrice * Math::SpecialFunctions::normalCDF( lDPlus ) -
             inStrike * Math::SpecialFunctions::normalCDF( lDMinus ) );
}

double funcBlackNegative( double inStrike, double inPrice, double inVol )
{
    if ( inVol == 0.0 ) { return std::max( { inStrike - inPrice, 0.0 } ); }
    double lDPlus  = std::log( inPrice / inStrike ) / inVol + 0.5 * inVol;
    double lDMinus = lDPlus - inVol;
    return ( -inPrice * Math::SpecialFunctions::normalCDF( -lDPlus ) +
             inStrike * Math::SpecialFunctions::normalCDF( -lDMinus ) );
}

std::vector<double> Model::calcZCBFromForwardRate(
    const std::vector<double>& inForwardRate ) const
{
    std::vector<double> lZCB( inForwardRate.size() );
    for ( std::size_t i = 0; i < lZCB.size(); ++i )
    {
        lZCB[i] =
            1.0 / ( 1.0 + ( mTerms[i + 1] - mTerms[i] ) * inForwardRate[i] );
    }
    for ( std::size_t i = 1; i < lZCB.size(); ++i ) { lZCB[i] *= lZCB[i - 1]; }
    return lZCB;
}
std::vector<double> Model::calcForwardRateFromZCB(
    const std::vector<double>& inZCB ) const
{
    std::vector<double> lForwardRate( inZCB.size() );
    lForwardRate[0] =
        ( 1.0 - inZCB[0] ) / ( inZCB[0] * ( mTerms[1] - mTerms[0] ) );
    for ( std::size_t i = 1; i < lForwardRate.size(); ++i )
    {
        lForwardRate[i] = ( inZCB[i - 1] - inZCB[i] ) /
                          ( inZCB[i] * ( mTerms[i + 1] - mTerms[i] ) );
    }
    return lForwardRate;
}
std::vector<double> Model::calcSumTauZCBFromZCB(
    const std::vector<double>& inZCB ) const
{
    std::vector<double> lSumTauZCB( inZCB.size() + 1, 0.0 );
    for ( std::size_t i = 0; i < inZCB.size(); ++i )
    {
        lSumTauZCB[i + 1] =
            lSumTauZCB[i] + inZCB[i] * ( mTerms[i + 1] - mTerms[i] );
    }
    return lSumTauZCB;
}

Model::Model( const std::vector<double>& inTerms ) : mTerms( inTerms )
{
    if ( mTerms.at( 0 ) != 0.0 )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::Model()\n" ) +
            std::string( "Terms.at(0) must be 0.0." ) );
    }
    if ( mTerms.size() <= 1 )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::Model()\n" ) +
            std::string( "Terms must have more than two components." ) );
    }
}
void Model::setInitForwardRate( const std::vector<double>& inForwardRate )
{
    if ( inForwardRate.size() != mTerms.size() - 1 )
    {
        throw std::invalid_argument(
            std::string(
                "Analytical::Black76::Model::setInitForwardRate()\n" ) +
            std::string( "size of forward rate must be Terms.size()-1." ) );
    }
    mForwardRate = inForwardRate;
    mZCB         = calcZCBFromForwardRate( inForwardRate );
    mSumTauZCB   = calcSumTauZCBFromZCB( mZCB );
}
void Model::setInitZCB( const std::vector<double>& inZCB )
{
    if ( inZCB.size() != mTerms.size() - 1 )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::setInitZCB()\n" ) +
            std::string( "size of ZCB must be Terms.size()-1." ) );
    }
    mZCB         = inZCB;
    mForwardRate = calcForwardRateFromZCB( inZCB );
    mSumTauZCB   = calcSumTauZCBFromZCB( mZCB );
}
void Model::setVol( double inVol ) { mVol = inVol; }
double Model::priceCaplet( double inStrike, std::size_t inIndex ) const
{
    if ( inIndex > mTerms.size() - 1 )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::priceCaplet()\n" ) +
            std::string( "invalid index of terms." ) );
    }
    return mZCB[inIndex] * ( mTerms[inIndex + 1] - mTerms[inIndex] ) *
           funcBlackPositive( inStrike, mForwardRate[inIndex],
                              mVol * std::sqrt( mTerms[inIndex] ) );
}
double Model::priceFloorlet( double inStrike, std::size_t inIndex ) const
{
    if ( inIndex > mTerms.size() - 1 )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::priceFloorlet()\n" ) +
            std::string( "invalid index of terms." ) );
    }
    return mZCB[inIndex] * ( mTerms[inIndex + 1] - mTerms[inIndex] ) *
           funcBlackNegative( inStrike, mForwardRate[inIndex],
                              mVol * std::sqrt( mTerms[inIndex] ) );
}
double Model::priceCap( double inStrike, std::size_t inIndStart,
                        std::size_t inIndLast ) const
{
    if ( inIndStart > inIndLast || inIndLast >= mTerms.size() )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::priceCap()\n" ) +
            std::string( "invalid index of terms." ) );
    }
    if ( inIndStart == inIndLast ) { return 0.0; }
    double lResult = 0.0;
    for ( std::size_t i = inIndStart; i < inIndLast; ++i )
    {
        lResult += priceCaplet( inStrike, i );
    }
    return lResult;
}
double Model::priceFloor( double inStrike, std::size_t inIndStart,
                          std::size_t inIndLast ) const
{
    if ( inIndStart > inIndLast || inIndLast >= mTerms.size() )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::priceFloor()\n" ) +
            std::string( "invalid index of terms." ) );
    }
    if ( inIndStart == inIndLast ) { return 0.0; }
    double lResult = 0.0;
    for ( std::size_t i = inIndStart; i < inIndLast; ++i )
    {
        lResult += priceFloorlet( inStrike, i );
    }
    return lResult;
}
double Model::priceSwapRate( std::size_t inIndStart,
                             std::size_t inIndLast ) const
{
    if ( inIndStart >= inIndLast || inIndLast >= mTerms.size() )
    {
        throw std::invalid_argument(
            std::string( "Analytical::Black76::Model::priceSwapRate()\n" ) +
            std::string( "invalid index of terms." ) );
    }
    double lZCBStart = ( inIndStart == 0 ) ? 1.0 : mZCB[inIndStart - 1];
    return ( lZCBStart - mZCB[inIndLast - 1] ) /
           ( mSumTauZCB[inIndLast] - mSumTauZCB[inIndStart] );
}
double Model::pricePayerSwaption( double inStrike, std::size_t inIndStart,
                                  std::size_t inIndLast ) const
{
    if ( inIndStart >= inIndLast || inIndLast >= mTerms.size() )
    {
        throw std::invalid_argument(
            std::string(
                "Analytical::Black76::Model::pricePayerSwaption()\n" ) +
            std::string( "invalid index of terms." ) );
    }
    return funcBlackPositive( inStrike, priceSwapRate( inIndStart, inIndLast ),
                              mVol * sqrt( mTerms[inIndStart] ) ) *
           ( mSumTauZCB[inIndLast] - mSumTauZCB[inIndStart] );
}
double Model::priceReceiverSwaption( double inStrike, std::size_t inIndStart,
                                     std::size_t inIndLast ) const
{
    if ( inIndStart >= inIndLast || inIndLast >= mTerms.size() )
    {
        throw std::invalid_argument(
            std::string(
                "Analytical::Black76::Model::priceReceiverSwaption()\n" ) +
            std::string( "invalid index of terms." ) );
    }
    return funcBlackNegative( inStrike, priceSwapRate( inIndStart, inIndLast ),
                              mVol * sqrt( mTerms[inIndStart] ) ) *
           ( mSumTauZCB[inIndLast] - mSumTauZCB[inIndStart] );
}

}  // namespace Black76
}  // namespace Analytical
