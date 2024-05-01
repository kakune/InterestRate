/**
 * @file vol_generator.cpp
 * @brief This implements volatility generator classes for forward model.
 * @author kakune
 * @date 4/15/2024
 */

#ifdef NINCLUDE_TPP
#include "LIBOR/forward/vol_generator.hpp"
#endif

#include <limits>

namespace LIBOR::Forward::VolGen
{

inline Constant::Constant( const Math::Vec& inVol ) : mVol( inVol ) {}

inline Math::Vec Constant::operator()(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inForwardRates,
    const Math::Vec& inStdBrownForFR ) const
{
    return mVol;
}

template <Process::RandomVec::C_StdBrown StdBrownVecGenerator_>
SABR<StdBrownVecGenerator_>::SABR( const Math::Vec& inInitVol,
                                   double inExponent, double inVolVol,
                                   const Math::Vec& inCorrSV,
                                   std::size_t inNPath,
                                   const Process::MarketData::Terms& inTerms,
                                   const Process::MarketData::Tenor& inTenor ) :
    msVols( std::make_shared<std::vector<Math::Vec>>( inNPath, inInitVol ) ),
    mExponent( inExponent ),
    mVolVol( inVolVol ),
    mCorrSV( inCorrSV ),
    mCorrSVInv( sqrt( 1.0 - inCorrSV * inCorrSV ) ),
    mTmpIndTerm( std::numeric_limits<std::size_t>::quiet_NaN() ),
    mVolVolSqrtDt( std::numeric_limits<double>::quiet_NaN() ),
    mTerms( inTerms ),
    mTenor( inTenor ),
    msStdBrownGen( std::make_shared<StdBrownVecGenerator_>( inInitVol.size() ) )
{
}

template <Process::RandomVec::C_StdBrown StdBrownVecGenerator_>
Math::Vec SABR<StdBrownVecGenerator_>::operator()(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inForwardRates,
    const Math::Vec& inStdBrownForFR ) const
{
    if ( mTmpIndTerm != inIndTerm )
    {
        msStdBrownGen->initialize( mTenor.minIndex( inIndTerm ) );
        mTmpIndTerm   = inIndTerm;
        mVolVolSqrtDt = mVolVol * mTerms.sqrtDifTime( inIndTerm );
    }
    Math::Vec lResult =
        ( *msVols )[inIndPath] *
        pow( inForwardRates[inIndPath][inIndTerm - 1], mExponent );
    Math::Vec lBrownVol =
        mCorrSV * inStdBrownForFR + mCorrSVInv * ( *msStdBrownGen )();
    ( *msVols )[inIndPath] *= ( 1.0 + mVolVolSqrtDt * lBrownVol );
    return lResult;
}

}  // namespace LIBOR::Forward::VolGen