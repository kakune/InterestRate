/**
 * @file Black_Scholes.cpp
 * @brief This implements analytical calculators relating to Black Scholes
 * Model.
 * @author kakune
 * @date 1/29/2024
 */

#include "analytical/Black_Scholes.hpp"

#include <cmath>

#include "math/special_functions.hpp"

namespace Analytical
{
namespace BlackScholes
{

double Model::priceEuropeanCallOption()
{
    double lD1 = ( std::log( mInitS / mStrike ) +
                   ( mRate + 0.5 * mVol * mVol ) * mTMat ) /
                 ( mVol * std::sqrt( mTMat ) );
    double lD2 = lD1 - mVol * std::sqrt( mTMat );
    return mInitS * Math::SpecialFunctions::normalCDF( lD1 ) -
           mStrike * std::exp( -mRate * mTMat ) *
               Math::SpecialFunctions::normalCDF( lD2 );
}
double Model::priceEuropeanPutOption()
{
    double lD1 = ( std::log( mInitS / mStrike ) +
                   ( mRate + 0.5 * mVol * mVol ) * mTMat ) /
                 ( mVol * std::sqrt( mTMat ) );
    double lD2 = lD1 - mVol * std::sqrt( mTMat );
    return -mInitS * Math::SpecialFunctions::normalCDF( -lD1 ) +
           mStrike * std::exp( -mRate * mTMat ) *
               Math::SpecialFunctions::normalCDF( -lD2 );
}

double priceEuropeanCallOption( double inInitS, double inStrike, double inRate,
                                double inVol, double inTMat )
{
    Model lObj{ inInitS, inStrike, inRate, inVol, inTMat };
    return lObj.priceEuropeanCallOption();
}
double priceEuropeanPutOption( double inInitS, double inStrike, double inRate,
                               double inVol, double inTMat )
{
    Model lObj{ inInitS, inStrike, inRate, inVol, inTMat };
    return lObj.priceEuropeanPutOption();
}

}  // namespace BlackScholes
}  // namespace Analytical
