/**
 * @file LFM.cpp
 * @brief This implements class for lognormal forward model.
 * @author kakune
 * @date 3/29/2024
 */

#include "LIBOR/forward/LFM.hpp"

namespace LIBOR
{
namespace Forward
{

Math::Vec ConstantLFM::driftTerm(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inStates ) const
{
    Math::Vec lDriftVec( mNFR );
    for ( std::size_t j = 0; j < mNFR; ++j )
    {
        double lFactor =
            std::exp( inStates[inIndPath][inIndTerm - 1]( j ) ) * mTauTenor[j];
        lDriftVec( j ) = lFactor / ( 1.0 + lFactor );
    }
    lDriftVec *= mVols;
    Math::Vec lDriftFactor( mNFR, 0.0 );
    for ( std::size_t k = 0; k < mNFR; ++k )
    {
        if ( mIndTenor[k + 1] < inIndTerm ) { continue; }
        for ( std::size_t j = 0; j <= k; ++j )
        {
            lDriftFactor( k ) += mRho( k, j ) * lDriftVec( j );
        }
    }
    return ( lDriftFactor * mVols + mNegativeHalfVolVol ) *
           mTerms.difTime( inIndTerm );
}
Math::Vec ConstantLFM::volTerm(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inStates,
    const Math::Vec& inRandomVec ) const
{
    return mVols * dot( mCorr, inRandomVec );
}
Math::Vec ConstantLFM::transfFRToState( const Math::Vec& inFR,
                                        std::size_t inIndTime ) const
{
    return log( inFR );
}
Math::Vec ConstantLFM::transfStateToFR( const Math::Vec& inState,
                                        std::size_t inIndTime ) const
{
    return exp( inState );
}

}  // namespace Forward
}  // namespace LIBOR
