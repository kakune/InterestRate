/**
 * @file integral_1d.tpp
 * @brief This implements the 1d-integral functions
 * @author kakune
 * @date 5/3/2024
 */

#define MATH_INTEGRAL_1D_MAX_ITER 1e8
#include <algorithm>
#include <cmath>
#ifdef NINCLUDE_TPP
#include "math/integral_1d.hpp"
#endif

namespace Math::Integral
{

namespace FiniteInterval
{

template <typename TypeX_, typename TypeY_>
TypeY_ trapz( auto inFunc, TypeX_ inMin, TypeX_ inMax, std::size_t inNDivision )
{
    const TypeX_ lDx = ( inMax - inMin ) / TypeX_( inNDivision );
    TypeX_ lTmpX     = inMin;
    TypeY_ lResult   = 0.5 * ( inFunc( inMin ) + inFunc( inMax ) );
    for ( std::size_t iX = 1; iX < inNDivision; ++iX )
    {
        lTmpX += lDx;
        lResult += inFunc( lTmpX );
    }
    lResult *= lDx;
    return lResult;
}

}  // namespace FiniteInterval

// https://www.st.nanzan-u.ac.jp/info/gr-thesis/ms/2010/05mm044.pdf
static inline double weightForOura( double inX )
{
    const double lCosh = cosh( M_PI_2 * sinh( inX ) );
    return 2.0 / ( lCosh * lCosh );
}
template <typename TypeY_>
static TypeY_ DEFormulaOura( auto inFunc, auto inTransfFunc, double inMin,
                             double inMax, double inTolRel,
                             std::size_t inInitNDivision = 25,
                             double inSafetyFactor       = 0.002 )
{
    TypeY_ lIntegralI, lIntegralJ, lFineIntegralI, lFineIntegralJ;
    // edge
    {
        const auto [lPhiAtMin, lDPhiAtMin] = inTransfFunc( inMin );
        const TypeY_ lFPhiAtMin            = inFunc( lPhiAtMin );
        lIntegralI                         = lFPhiAtMin * lDPhiAtMin;
        lIntegralJ = lFPhiAtMin * weightForOura( inMin );

        const auto [lPhiAtMax, lDPhiAtMax] = inTransfFunc( inMax );
        const TypeY_ lFPhiAtMax            = inFunc( lPhiAtMax );
        lIntegralI += lFPhiAtMax * lDPhiAtMax;
        lIntegralJ += lFPhiAtMax * weightForOura( inMax );

        lIntegralI *= 0.5;
        lIntegralJ *= 0.5;
    }

    std::size_t lTmpNDivision = inInitNDivision;
    double lDifX              = ( inMax - inMin ) / double( lTmpNDivision );
    {
        double lTmpX = inMin;
        for ( std::size_t iDiv = 1; iDiv < lTmpNDivision; ++iDiv )
        {
            lTmpX += lDifX;
            const auto [lPhi, lDPhi] = inTransfFunc( lTmpX );
            const TypeY_ lFPhi       = inFunc( lPhi );
            lIntegralI += lFPhi * lDPhi;
            lIntegralJ += lFPhi * weightForOura( lTmpX );
        }
    }
    lIntegralI *= lDifX;
    lIntegralJ *= lDifX;

    const double inBound = inSafetyFactor * sqrt( inTolRel );
    do {
        double lHalfDifX                         = 0.5 * lDifX;
        double lTmpX                             = inMin + lHalfDifX;
        const auto [lPhiAtTmpMin, lDPhiAtTmpMin] = inTransfFunc( lTmpX );
        const TypeY_ lFPhiAtTmpMin               = inFunc( lPhiAtTmpMin );
        TypeY_ lAdditionalI = lFPhiAtTmpMin * lDPhiAtTmpMin;
        TypeY_ lAdditionalJ = lFPhiAtTmpMin * weightForOura( lTmpX );
        for ( std::size_t iDiv = 1; iDiv < lTmpNDivision; ++iDiv )
        {
            lTmpX += lDifX;
            const auto [lPhi, lDPhi] = inTransfFunc( lTmpX );
            const TypeY_ lFPhi       = inFunc( lPhi );
            lAdditionalI += lFPhi * lDPhi;
            lAdditionalJ += lFPhi * weightForOura( lTmpX );
        }
        lFineIntegralI = 0.5 * lIntegralI + lHalfDifX * lAdditionalI;
        lFineIntegralJ = 0.5 * lIntegralJ + lHalfDifX * lAdditionalJ;
        std::swap( lIntegralI, lFineIntegralI );
        std::swap( lIntegralJ, lFineIntegralJ );
        std::swap( lDifX, lHalfDifX );
        lTmpNDivision += lTmpNDivision;
    } while ( lTmpNDivision < MATH_INTEGRAL_1D_MAX_ITER &&
              std::max( { std::abs( lFineIntegralI - lIntegralI ),
                          std::abs( lFineIntegralJ - lIntegralJ ) } ) >=
                  inBound * std::abs( lIntegralI ) );

    return lIntegralI;
}

namespace FiniteInterval
{

template <typename TypeY_>
TypeY_ DEFormula( auto inFunc, double inMin, double inMax, double inTolRel )
{
    const double lHalfInterval    = 0.5 * ( inMax - inMin );
    const double lHalfSumEndPoint = 0.5 * ( inMax + inMin );
    //! make interval -1~1
    auto lFuncForAdaptive = [&inFunc, lHalfInterval,
                             lHalfSumEndPoint]( double inX ) -> TypeY_
    {
        const double lTmpX = inX * lHalfInterval + lHalfSumEndPoint;
        return lHalfInterval * inFunc( lTmpX );
    };
    auto lTransfFunc = []( double inX ) -> std::pair<double, double>
    {
        const double lSinh = sinh( inX );
        const double lPhi  = tanh( M_PI_2 * lSinh );
        const double lDPhi =
            M_PI * cosh( inX ) / ( 1.0 + cosh( M_PI * lSinh ) );
        return { lPhi, lDPhi };
    };
    return DEFormulaOura<TypeY_>( lFuncForAdaptive, lTransfFunc, -3.1905,
                                  3.1905, inTolRel );
}

}  // namespace FiniteInterval

namespace InfiniteInterval
{

template <typename TypeY_> TypeY_ DEFormula( auto inFunc, double inTolRel )
{
    auto lTransfFunc = []( double inX ) -> std::pair<double, double>
    {
        const double lHalfPiSinh = M_PI_2 * sinh( inX );
        const double lPhi        = sinh( lHalfPiSinh );
        const double lDPhi       = M_PI_2 * cosh( inX ) * cosh( lHalfPiSinh );
        return { lPhi, lDPhi };
    };
    return DEFormulaOura<TypeY_>( inFunc, lTransfFunc, -6.0, 6.0, inTolRel );
}

}  // namespace InfiniteInterval

namespace UpperInfiniteInterval
{

template <typename TypeY_>
TypeY_ DEFormula( auto inFunc, double inMin, double inTolRel )
{
    // make interval 0~infty
    auto lFuncForAdaptive = [&inFunc, inMin]( double inX ) -> TypeY_
    { return inFunc( inX + inMin ); };
    auto lTransfFunc = []( double inX ) -> std::pair<double, double>
    {
        const double lPhi = exp( M_PI_2 * sinh( inX ) );
        return { lPhi, M_PI_2 * lPhi * cosh( inX ) };
    };
    return DEFormulaOura<TypeY_>( lFuncForAdaptive, lTransfFunc, -6.75, 6.75,
                                  inTolRel );
}

template <typename TypeY_>
TypeY_ DEFormulaForExp( auto inFunc, double inMin, double inTolRel )
{
    // make interval 0~infty
    auto lFuncForAdaptive = [&inFunc, inMin]( double inX ) -> TypeY_
    { return inFunc( inX + inMin ); };
    auto lTransfFunc = []( double inX ) -> std::pair<double, double>
    {
        const double lExp = exp( -inX );
        const double lPhi = exp( inX - lExp );
        return { lPhi, ( 1.0 + lExp ) * lPhi };
    };
    return DEFormulaOura<TypeY_>( lFuncForAdaptive, lTransfFunc, -5.0, 100.0,
                                  inTolRel );
}

template <typename TypeY_>
TypeY_ DEFormulaForGaussian( auto inFunc, double inMin, double inTolRel )
{
    // make interval 0~infty
    auto lFuncForAdaptive = [&inFunc, inMin]( double inX ) -> TypeY_
    { return inFunc( inX + inMin ); };
    auto lTransfFunc = []( double inX ) -> std::pair<double, double>
    {
        const double lExp = std::exp( -inX );
        const double lPhi = std::exp( 0.5 * inX - lExp );
        return { lPhi, lPhi * ( 0.5 + lExp ) };
    };
    return DEFormulaOura<TypeY_>( lFuncForAdaptive, lTransfFunc, -6.0, 100.0,
                                  inTolRel );
}

template <typename TypeY_>
TypeY_ DEFormulaForOscillation( auto inFunc, double inMin, double inTolRel,
                                double inFactorK )
{
    // make interval 0~infty
    auto lFuncForAdaptive = [&inFunc, inMin]( double inX ) -> TypeY_
    { return inFunc( inX + inMin ); };
    auto lTransfFunc = [inFactorK]( double inX ) -> std::pair<double, double>
    {
        if ( inX > 5.0 ) { return { inX, 1.0 }; }
        const double lExpSinh = exp( -inFactorK * sinh( inX ) );
        const double lDenom   = 1.0 / ( 1.0 - lExpSinh );
        return { inX * lDenom,
                 ( 1.0 - ( 1.0 + inFactorK * inX * cosh( inX ) ) * lExpSinh ) *
                     lDenom * lDenom };
    };
    return DEFormulaOura<TypeY_>( lFuncForAdaptive, lTransfFunc, -5.0, 1e6,
                                  inTolRel );
}

}  // namespace UpperInfiniteInterval

}  // namespace Math::Integral
