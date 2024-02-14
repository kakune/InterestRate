/**
 * @file asset.cpp
 * @brief This implements calculator using CUDA of assets.
 * @author kakune
 * @date 2/14/2024
 */

#include "process_cuda/asset.hpp"

#include <memory>
#include <vector>

#include "analytical/Black_Scholes.hpp"
#include "math/findroot_1d.hpp"

namespace ProcessCUDA
{
namespace Asset
{

double ModelForwardAbstract::impliedVolatility( double inStrike )
{
    // if ( inStrike > mInitPrice )
    {
        double lCallPrice = priceCallOption( inStrike );
        auto lFunction = [this, inStrike, lCallPrice]( double inVol ) -> double
        {
            return lCallPrice -
                   Analytical::BlackScholes::europeanCallOptionPrice(
                       this->mInitPrice, inStrike, 0.0, inVol,
                       this->mTimeMaturity );
        };
        return Math::FindRoot1D::Brent( lFunction, 0.000000001, 10.0 );
    }
    double lPutPrice = pricePutOption( inStrike );
    auto lFunction   = [this, inStrike, lPutPrice]( double inVol ) -> double
    {
        return lPutPrice - Analytical::BlackScholes::europeanPutOptionPrice(
                               this->mInitPrice, inStrike, 0.0, inVol,
                               this->mTimeMaturity );
    };
    return Math::FindRoot1D::Brent( lFunction, 0.000000001, 10.0 );
}

}  // namespace Asset
}  // namespace ProcessCUDA