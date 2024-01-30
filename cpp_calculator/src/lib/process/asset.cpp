/**
 * @file asset.hpp
 * @brief This implements calculator of assets.
 * @author kakune
 * @date 1/30/2024
 */

#include "process/asset.hpp"

#include <memory>
#include <vector>

#include "analytical/Black_Scholes.hpp"
#include "math/findroot_1d.hpp"

namespace Process
{
namespace Asset
{

double ModelForwardAbstract::priceCallOption( double inStrike,
                                              std::size_t inIndTime )
{
    double lResult = 0;
    for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
    {
        lResult += std::max(
            0.0, mForwardPrice.at( iPath ).at( inIndTime ) - inStrike );
    }
    return ( lResult / double( mNPath ) );
}

double ModelForwardAbstract::impliedVolatility( double inStrike,
                                                std::size_t inIndTime )
{
    double lCallPrice    = priceCallOption( inStrike, inIndTime );
    double lTimeMaturity = msTerms->at( inIndTime ) - msTerms->front();
    auto lFunction       = [this, inStrike, lCallPrice,
                      lTimeMaturity]( double inVol ) -> double
    {
        return lCallPrice -
               Analytical::BlackScholes::europeanCallOptionPrice(
                   this->mInitPrice, inStrike, 0.0, inVol, lTimeMaturity );
    };
    return Math::FindRoot1D::Brent( lFunction, 0.000001, 10.0 );
}

void LocalVolatilityForwardAbstract::calcEachForwardPrice()
{
    muRandomPath->makeRandomVariables();
    for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
    {
        for ( std::size_t iTime = 1; iTime < msTerms->size(); ++iTime )
        {
            mForwardPrice.at( iPath ).at( iTime ) =
                mForwardPrice.at( iPath ).at( iTime - 1 ) +
                muRandomPath->at( iPath ).at( iTime ) *
                    localVolatility( mForwardPrice.at( iPath ).at( iTime - 1 ),
                                     msTerms->at( iTime - 1 ) );
        }
    }
}

double BlackSholesForward::localVolatility( double inPrice, double inTime )
{
    return mVol * std::pow( inPrice, 0.8 );
}

}  // namespace Asset
}  // namespace Process