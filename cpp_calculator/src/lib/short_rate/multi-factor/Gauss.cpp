/**
 * @file Gauss.cpp
 * @brief This implements multi-factor Gaussian short-rate models. (see Chap.4
 * of Brigo)
 * @author kakune
 * @date 3/24/2024
 */

#include "short_rate/multi-factor/Gauss.hpp"

#include <iostream>
#include <memory>

namespace ShortRate
{
namespace MultiFactor
{

Math::Vec ConstantGauss::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots ) const
{
    return mDriftCoeff * inSpots[inIndPath][inIndTerm - 1];
}
Math::Vec ConstantGauss::volTerm(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots,
    const Math::Vec& inRandomVec ) const
{
    return dotLowerMatVec( mVolCoeff, inRandomVec );
}
double ConstantGauss::transfStateToRate( const Math::Vec& inState,
                                         std::size_t inIndTime ) const
{
    return inState.sum();
}

double ConstantGauss::analyticalPriceZCB( double inTerminalTime ) const
{
    if ( mDim != 2 )
    {
        std::cerr
            << "Error : "
               "Process::ShortRateMCMulti::ConstantGauss::analyticalPriceZCB()"
            << std::endl
            << "Sorry, Dim != 2 is not implemented....." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    double lTime     = inTerminalTime - mTerms[0];
    double lA        = -mDriftCoeff( 0 );
    double lExpA     = std::exp( -lA * lTime );
    double lB        = -mDriftCoeff( 1 );
    double lExpB     = std::exp( -lB * lTime );
    double lSumAB    = lA + lB;
    double lAB       = lA * lB;
    double lSigma    = mVolCoeff( 0, 0 );
    double lEta      = std::sqrt( mVolCoeff( 1, 0 ) * mVolCoeff( 1, 0 ) +
                                  mVolCoeff( 1, 1 ) * mVolCoeff( 1, 1 ) );
    double lRho      = mVolCoeff( 1, 0 ) / lEta;
    double lFactorA  = ( 1.0 - lExpA ) / lA;
    double lFactorB  = ( 1.0 - lExpB ) / lB;
    double lFactorAB = ( 1.0 - std::exp( -lSumAB * lTime ) ) / lSumAB;
    double lV = 2.0 * ( lTime - lFactorA - lFactorB + lFactorAB ) * lRho *
                lSigma * lEta / lAB;
    lV += ( lSigma * lSigma / ( lA * lA ) ) *
          ( lTime + ( 2.0 * lExpA - 0.5 * lExpA * lExpA - 1.5 ) / lA );
    lV += ( lEta * lEta / ( lB * lB ) ) *
          ( lTime + ( 2.0 * lExpB - 0.5 * lExpB * lExpB - 1.5 ) / lB );
    return std::exp( -lFactorA * mInitState( 0 ) - lFactorB * mInitState( 1 ) +
                     0.5 * lV );
}

Math::Vec G2ppWithMarket::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots ) const
{
    return mDriftCoeff * inSpots[inIndPath][inIndTerm - 1];
}
Math::Vec G2ppWithMarket::volTerm(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots,
    const Math::Vec& inRandomVec ) const
{
    return dotLowerMatVec( mVolCoeff, inRandomVec );
}
double G2ppWithMarket::transfStateToRate( const Math::Vec& inState,
                                          std::size_t inIndTime ) const
{
    double lDifTime = mTerms[inIndTime] - mTerms[0];
    double lExpA    = 1.0 - std::exp( -mA * lDifTime );
    double lExpB    = 1.0 - std::exp( -mB * lDifTime );
    return inState.sum() + mInstantaneousFRs[inIndTime] +
           mFactorA * lExpA * lExpA + mFactorB * lExpB * lExpB +
           mFactorAB * lExpA * lExpB;
    // return inState.sum() +
    //        mMarketZCB.instantaneousForwardRate( mTerms[inIndTime] ) +
    //        mFactorA * lExpA * lExpA + mFactorB * lExpB * lExpB +
    //        mFactorAB * lExpA * lExpB;
}
void G2ppWithMarket::prepareFactors()
{
    mA            = -mDriftCoeff( 0 );
    mB            = -mDriftCoeff( 1 );
    double lSigma = mVolCoeff( 0, 0 );
    double lEta   = std::sqrt( mVolCoeff( 1, 0 ) * mVolCoeff( 1, 0 ) +
                               mVolCoeff( 1, 1 ) * mVolCoeff( 1, 1 ) );
    double lRho   = mVolCoeff( 1, 0 ) / lEta;
    mFactorA      = 0.5 * lSigma * lSigma / ( mA * mA );
    mFactorB      = 0.5 * lEta * lEta / ( mB * mB );
    mFactorAB     = lRho * lSigma * lEta / ( mA * mB );
    for ( std::size_t i = 0; i < mTerms.size(); ++i )
    {
        mInstantaneousFRs[i] = mMarketZCB.instantaneousForwardRate( mTerms[i] );
    }
}

}  // namespace MultiFactor
}  // namespace ShortRate
