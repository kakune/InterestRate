/**
 * @file one-factor_Affine.cpp
 * @brief This implements one-factor Affine short-rate models. (see 3.2.4 of
 * Brigo)
 * @author kakune
 * @date 3/14/2024
 */

#include "process/short_rate/one-factor_Affine.hpp"

#include <iostream>
#include <memory>
#include <vector>

#include "process/market.hpp"

namespace Process
{
namespace ShortRate
{

double ConstantAffine::driftCoeff( std::size_t inIndPath,
                                   std::size_t inIndTerm ) const
{
    return mLambda * mSpotRates.at( inIndPath ).at( inIndTerm - 1 ) + mEta;
}
double ConstantAffine::volCoeff( std::size_t inIndPath,
                                 std::size_t inIndTerm ) const
{
    return std::sqrt( mGamma * mSpotRates.at( inIndPath ).at( inIndTerm - 1 ) +
                      mDelta );
}

}  // namespace ShortRate
}  // namespace Process
