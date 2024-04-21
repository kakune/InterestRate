/**
 * @file step_calculator.tpp
 * @brief This implements step calculator classes for forward model.
 * @author kakune
 * @date 4/15/2024
 */

#ifndef LIBOR_FORWARD_STEP_CALCULATOR_TPP
#define LIBOR_FORWARD_STEP_CALCULATOR_TPP

#include "LIBOR/forward/step_calculator.hpp"

namespace LIBOR::Forward::StepCalc
{

Math::Mat calcUpperTriangularRhoWithOutDiag( const Math::Mat& inCorr )
{
    Math::Mat lResult = dot( inCorr, inCorr.transpose() );
    for ( std::size_t i = 0; i < lResult.sizeRow(); ++i )
    {
        for ( std::size_t j = 0; j <= i; ++j ) { lResult( i, j ) = 0.0; }
    }
    return lResult;
}

template <C_VolGen VolatilityGenerator_>
LogNormalTerminalMeas<VolatilityGenerator_>::LogNormalTerminalMeas(
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::Tenor& inTenor, const Math::Mat& inCorr,
    const VolatilityGenerator_& inVol ) :
    mTerms( inTerms ),
    mTenor( inTenor ),
    mCorr( inCorr ),
    mRho( dot( inCorr, inCorr.transpose() ) ),
    mUpperRho( calcUpperTriangularRhoWithOutDiag( inCorr ) ),
    mVolGen( inVol )
{
}

template <C_VolGen VolatilityGenerator_>
Math::Vec LogNormalTerminalMeas<VolatilityGenerator_>::operator()(
    std::size_t inIndPath, std::size_t inIndTerm, const StatesType& inStates,
    const Math::Vec& inStdBrownianVec ) const
{
    if ( mTenor.minIndex( inIndTerm ) >= mTenor.size() )
    {
        return Math::Vec( mTenor.size(), 0.0 );
    }
    Math::Vec lVol =
        mVolGen( inIndPath, inIndTerm, inStates.getValues(), inStdBrownianVec );
    Math::Vec lTauFR =
        mTenor.getTauVec() * inStates.getValues()[inIndPath][inIndTerm - 1];
    Math::Vec lDriftFactor = lVol * lTauFR / ( 1.0 + lTauFR );
    Math::Vec lDrift       = -dot( mUpperRho, lDriftFactor ) * lVol *
                       inStates.getValues()[inIndPath][inIndTerm - 1] *
                       mTerms.difTime( inIndTerm );
    Math::Vec lResult =
        lDrift + lVol * inStates.getValues()[inIndPath][inIndTerm - 1] *
                     mTerms.sqrtDifTime( inIndTerm ) *
                     dot( mCorr, inStdBrownianVec );
    for ( std::size_t i = 0; i < mTenor.minIndex( inIndTerm ); ++i )
    {
        lResult( i ) = 0.0;
    }
    return lResult;
}

template <C_VolGen VolatilityGenerator_>
LogNormalTerminalMeasWithLog<VolatilityGenerator_>::
    LogNormalTerminalMeasWithLog( const Process::MarketData::Terms& inTerms,
                                  const Process::MarketData::Tenor& inTenor,
                                  const Math::Mat& inCorr,
                                  const VolatilityGenerator_& inVol ) :
    mTerms( inTerms ),
    mTenor( inTenor ),
    mCorr( inCorr ),
    mRho( dot( inCorr, inCorr.transpose() ) ),
    mUpperRho( calcUpperTriangularRhoWithOutDiag( inCorr ) ),
    mVolGen( inVol )
{
}

template <C_VolGen VolatilityGenerator_>
Math::Vec LogNormalTerminalMeasWithLog<VolatilityGenerator_>::operator()(
    std::size_t inIndPath, std::size_t inIndTerm, const StatesType& inStates,
    const Math::Vec& inStdBrownianVec ) const
{
    if ( mTenor.minIndex( inIndTerm ) >= mTenor.size() )
    {
        return Math::Vec( mTenor.size(), 0.0 );
    }
    Math::Vec lVol =
        mVolGen( inIndPath, inIndTerm, inStates.getValues(), inStdBrownianVec );
    Math::Vec lTauFR =
        mTenor.getTauVec() * inStates.getValues()[inIndPath][inIndTerm - 1];
    Math::Vec lDriftFactor = lVol * lTauFR / ( 1.0 + lTauFR );
    Math::Vec lDrift = ( -dot( mUpperRho, lDriftFactor ) - 0.5 * lVol ) * lVol *
                       mTerms.difTime( inIndTerm );
    Math::Vec lResult = lDrift + lVol * mTerms.sqrtDifTime( inIndTerm ) *
                                     dot( mCorr, inStdBrownianVec );
    for ( std::size_t i = 0; i < mTenor.minIndex( inIndTerm ); ++i )
    {
        lResult( i ) = 0.0;
    }
    return lResult;
}

template <C_VolGen VolatilityGenerator_>
NormalTerminalMeas<VolatilityGenerator_>::NormalTerminalMeas(
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::Tenor& inTenor, const Math::Mat& inCorr,
    const VolatilityGenerator_& inVol ) :
    mTerms( inTerms ),
    mTenor( inTenor ),
    mCorr( inCorr ),
    mRho( dot( inCorr, inCorr.transpose() ) ),
    mUpperRho( calcUpperTriangularRhoWithOutDiag( inCorr ) ),
    mVolGen( inVol )
{
}

template <C_VolGen VolatilityGenerator_>
Math::Vec NormalTerminalMeas<VolatilityGenerator_>::operator()(
    std::size_t inIndPath, std::size_t inIndTerm, const StatesType& inStates,
    const Math::Vec& inStdBrownianVec ) const
{
    if ( mTenor.minIndex( inIndTerm ) >= mTenor.size() )
    {
        return Math::Vec( mTenor.size(), 0.0 );
    }
    Math::Vec lVol =
        mVolGen( inIndPath, inIndTerm, inStates.getValues(), inStdBrownianVec );

    Math::Vec lDriftFactor =
        lVol * mTenor.getTauVec() /
        ( 1.0 +
          mTenor.getTauVec() * inStates.getValues()[inIndPath][inIndTerm - 1] );
    Math::Vec lDrift =
        -dot( mUpperRho, lDriftFactor ) * lVol * mTerms.difTime( inIndTerm );
    Math::Vec lResult = lDrift + lVol * mTerms.sqrtDifTime( inIndTerm ) *
                                     dot( mCorr, inStdBrownianVec );
    for ( std::size_t i = 0; i < mTenor.minIndex( inIndTerm ); ++i )
    {
        lResult( i ) = 0.0;
    }
    return lResult;
}

}  // namespace LIBOR::Forward::StepCalc

#endif