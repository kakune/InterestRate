/**
 * @file step_calculator.tpp
 * @brief This implements step calculator classes for forward model.
 * @author kakune
 * @date 4/15/2024
 */

#ifdef NINCLUDE_TPP
#include "LIBOR/forward/model_data.tpp"
#include "LIBOR/forward/step_calculator.hpp"
#endif

namespace LIBOR::Forward::StepCalc
{

static inline Math::Mat calcUpperTriangularRhoWithOutDiag(
    const Math::Mat& inCorr )
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
    Math::Vec lVol = mVolGen( inIndPath, inIndTerm, inStates.getForwards(),
                              inStdBrownianVec );
    Math::Vec lTauFR =
        mTenor.getTauVec() * inStates.getForwards()[inIndPath][inIndTerm - 1];
    Math::Vec lDriftFactor = lVol * lTauFR / ( 1.0 + lTauFR );
    Math::Vec lDrift       = -dot( mUpperRho, lDriftFactor ) * lVol *
                       inStates.getForwards()[inIndPath][inIndTerm - 1] *
                       mTerms.difTime( inIndTerm );
    Math::Vec lResult =
        lDrift + lVol * inStates.getForwards()[inIndPath][inIndTerm - 1] *
                     mTerms.sqrtDifTime( inIndTerm ) *
                     dot( mCorr, inStdBrownianVec );
    for ( std::size_t i = 0; i < mTenor.minIndex( inIndTerm ); ++i )
    {
        lResult[i] = 0.0;
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
    Math::Vec lVol = mVolGen( inIndPath, inIndTerm, inStates.getForwards(),
                              inStdBrownianVec );
    Math::Vec lTauFR =
        mTenor.getTauVec() * inStates.getForwards()[inIndPath][inIndTerm - 1];
    Math::Vec lDriftFactor = lVol * lTauFR / ( 1.0 + lTauFR );
    Math::Vec lDrift = ( -dot( mUpperRho, lDriftFactor ) - 0.5 * lVol ) * lVol *
                       mTerms.difTime( inIndTerm );
    Math::Vec lResult = lDrift + lVol * mTerms.sqrtDifTime( inIndTerm ) *
                                     dot( mCorr, inStdBrownianVec );
    for ( std::size_t i = 0; i < mTenor.minIndex( inIndTerm ); ++i )
    {
        lResult[i] = 0.0;
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
    Math::Vec lVol = mVolGen( inIndPath, inIndTerm, inStates.getForwards(),
                              inStdBrownianVec );

    Math::Vec lDriftFactor =
        lVol * mTenor.getTauVec() /
        ( 1.0 + mTenor.getTauVec() *
                    inStates.getForwards()[inIndPath][inIndTerm - 1] );
    Math::Vec lDrift =
        -dot( mUpperRho, lDriftFactor ) * lVol * mTerms.difTime( inIndTerm );
    Math::Vec lResult = lDrift + lVol * mTerms.sqrtDifTime( inIndTerm ) *
                                     dot( mCorr, inStdBrownianVec );
    for ( std::size_t i = 0; i < mTenor.minIndex( inIndTerm ); ++i )
    {
        lResult[i] = 0.0;
    }
    return lResult;
}

static inline std::vector<Math::Mat> calcSpotRhos(
    const Math::Mat& inCorr, const Process::MarketData::Tenor& inTenor )
{
    Math::Mat lRho = dot( inCorr.transpose(), inCorr );
    std::vector<Math::Mat> lResult( inTenor.size() + 1, lRho );
    for ( std::size_t i = 0; i < lResult.size(); ++i )
    {
        for ( std::size_t j = 0; j < i; ++j )
        {
            for ( std::size_t k = 0; k < lRho.sizeRow(); ++k )
            {
                lResult[i]( k, j ) = 0.0;
            }
        }
    }
    return lResult;
}

template <C_VolGen VolatilityGenerator_>
LogNormalSpotMeas<VolatilityGenerator_>::LogNormalSpotMeas(
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::Tenor& inTenor, const Math::Mat& inCorr,
    const VolatilityGenerator_& inVol ) :
    mTerms( inTerms ),
    mTenor( inTenor ),
    mCorr( inCorr ),
    mRho( dot( inCorr, inCorr.transpose() ) ),
    mStepRhos( calcSpotRhos( inCorr, inTenor ) ),
    mVolGen( inVol )
{
}

template <C_VolGen VolatilityGenerator_>
Math::Vec LogNormalSpotMeas<VolatilityGenerator_>::operator()(
    std::size_t inIndPath, std::size_t inIndTerm, const StatesType& inStates,
    const Math::Vec& inStdBrownianVec ) const
{
    std::size_t lMinIndex = mTenor.minIndex( inIndTerm );
    const Math::Vec& lFR  = inStates.getForwards()[inIndPath][inIndTerm - 1];
    Math::Vec lVol   = mVolGen( inIndPath, inIndTerm, inStates.getForwards(),
                                inStdBrownianVec );
    Math::Vec lTauFR = mTenor.getTauVec() * lFR;
    Math::Vec lDriftFactor = lVol * lTauFR / ( 1.0 + lTauFR );
    Math::Vec lDrift = dot( mStepRhos[lMinIndex], lDriftFactor ) * lVol * lFR *
                       mTerms.difTime( inIndTerm );
    Math::Vec lResult = lDrift + lVol * lFR * mTerms.sqrtDifTime( inIndTerm ) *
                                     dot( mCorr, inStdBrownianVec );
    for ( std::size_t i = 0; i < lMinIndex; ++i ) { lResult[i] = 0.0; }
    return lResult;
}

}  // namespace LIBOR::Forward::StepCalc