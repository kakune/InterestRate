/**
 * @file one-factor_Gauss.cpp
 * @brief This implements one-factor Gaussian short-rate models. (see Chap.10 of
 * Andersen & Piterbarg)
 * @author kakune
 * @date 3/8/2024
 */

#include "process/short_rate/one-factor_Gauss.hpp"

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

}  // namespace ShortRate
}  // namespace Process
