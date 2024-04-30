/**
 * @file optimize.hpp
 * @brief This implements optimizers.
 * @author kakune
 * @date 4/30/2024
 */

#ifdef NINCLUDE_TPP
#include "math/optimize.hpp"
#endif

#include <cmath>

namespace Math::Optimize
{

static void calcNumericalGradient( C_ObjectiveFunction auto inObjectiveFunc,
                                   const Math::Vec& inInitialValue,
                                   Math::Vec& inGrad, double inDif = 1e-8 )
{
    for ( size_t i = 0; i < inInitialValue.size(); ++i )
    {
        Math::Vec lVecPlusDif  = inInitialValue;
        Math::Vec lVecMinusDif = inInitialValue;
        lVecPlusDif[i] += inDif;
        lVecMinusDif[i] -= inDif;
        inGrad[i] = ( inObjectiveFunc( lVecPlusDif ) -
                      inObjectiveFunc( lVecMinusDif ) );
    }
    inGrad *= ( 0.5 / inDif );
}

Math::Vec adam( C_ObjectiveFunction auto inObjectiveFunc,
                const Math::Vec& inVec, double inAlpha, double inBeta1,
                double inBeta2, double inEpsilon, std::size_t inMaxIter,
                double inEpsGradNorm )
{
    Math::Vec lVelocity = Math::makeVec( inVec.size(), 0.0 );
    Math::Vec lAccel    = Math::makeVec( inVec.size(), 1.0 );
    Math::Vec lGrad     = Math::makeVec( inVec.size() );
    Math::Vec lResult   = Math::makeVec( inVec );

    double lPowBeta1 = 1.0;
    double lPowBeta2 = 1.0;
    double lDifBeta1 = 1.0 - inBeta1;
    double lDifBeta2 = 1.0 - inBeta2;
    for ( std::size_t iIter = 0; iIter < inMaxIter; ++iIter )
    {
        lPowBeta1 *= inBeta1;
        lPowBeta2 *= inBeta2;

        calcNumericalGradient( inObjectiveFunc, lResult, lGrad );
        lAccel    = inBeta1 * lAccel + lDifBeta1 * lGrad;
        lVelocity = inBeta2 * lVelocity + lDifBeta2 * lGrad * lGrad;

        lResult -= inAlpha * ( lAccel / ( 1.0 - lPowBeta1 ) ) /
                   ( sqrt( lVelocity / ( 1.0 - lPowBeta2 ) ) + inEpsilon );

        if ( Math::dot( lGrad, lGrad ) < inEpsGradNorm ) { break; }
    }
    return lResult;
}

}  // namespace Math::Optimize