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
               msTerms->at( inIndPath - 1 ), 1 ) +
           mVol2 * ( msTerms->at( inIndPath - 1 ) - msTerms->at( 0 ) );
}
double HoLee::volCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const
{
    return mVol;
}

}  // namespace ShortRate
}  // namespace Process
