/**
 * @file Black76.cpp
 * @brief This implements analytical calculators relating to Black-76 Model.
 * Model.
 * @author kakune
 * @date 4/13/2024
 */

#include "analytical/Black76.hpp"

#include <cmath>

#include "math/special_functions.hpp"

namespace Analytical
{
namespace Black76
{

double Model::operator()( bool inIsCaplet )
{
    double lOmega = 1.0;
    if ( !inIsCaplet ) { lOmega = -1.0; }
    double lV      = mVol * std::sqrt( mTimeStart );
    double lDPlus  = std::log( mCurrentFR / mStrike ) / lV + 0.5 * lV;
    double lDMinus = lDPlus - lV;
    return lOmega * mZCB * mTau *
           ( mCurrentFR * Math::SpecialFunctions::normalCDF( lOmega * lDPlus ) -
             mStrike * Math::SpecialFunctions::normalCDF( lOmega * lDMinus ) );
}
double Model::priceCaplet() { return ( *this )( true ); }
double Model::priceFloorlet() { return ( *this )( false ); }

double priceCaplet( double inCurrentFR, double inStrike, double inTimeStart,
                    double inTau, double inVol, double inZCB )
{
    Model lObj{ inCurrentFR, inStrike, inTimeStart, inTau, inVol, inZCB };
    return lObj.priceCaplet();
}
double priceFloorlet( double inCurrentFR, double inStrike, double inTimeStart,
                      double inTau, double inVol, double inZCB )
{
    Model lObj{ inCurrentFR, inStrike, inTimeStart, inTau, inVol, inZCB };
    return lObj.priceFloorlet();
}

}  // namespace Black76
}  // namespace Analytical
