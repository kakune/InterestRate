/**
 * @file Affine.cpp
 * @brief This implements multi-factor Affine short-rate models. (see Chap.4
 * of Brigo)
 * @author kakune
 * @date 3/25/2024
 */

#include "short_rate/multi-factor/Affine.hpp"

#include "math/matrix.hpp"

namespace ShortRate
{
namespace MultiFactor
{

Math::Vec CIR2ppWithMarket::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots ) const
{
    return mConvSH * ( mMean - inSpots[inIndPath][inIndTerm - 1] );
}
Math::Vec CIR2ppWithMarket::volTerm(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<Math::Vec>>& inSpots,
    const Math::Vec& inRandomVec ) const
{
    return mVol * sqrt( inSpots[inIndPath][inIndTerm] ) * inRandomVec;
}
double CIR2ppWithMarket::transfStateToRate( const Math::Vec& inState,
                                            std::size_t inIndTime ) const
{
    return inState.sum() + mShifts[inIndTime];
}
std::vector<double> CIR2ppWithMarket::calcShifts() const
{
    std::vector<double> lResult( mTerms.size() );
    Math::Vec lFactorH( sqrt( mConvSH * mConvSH + 2.0 * mVol * mVol ) );
    Math::Vec lTwoH     = 2.0 * lFactorH;
    Math::Vec lTwoKMean = 2.0 * mConvSH * mMean;
    Math::Vec lKpH      = mConvSH + lFactorH;
    Math::Vec lFourHHX  = lTwoH * lTwoH * mInitState;
    for ( std::size_t iTime = 0; iTime < mTerms.size(); ++iTime )
    {
        Math::Vec lExp     = exp( ( mTerms[iTime] - mTerms[0] ) * lFactorH );
        Math::Vec lExpmOne = lExp - 1.0;
        Math::Vec lDenom   = 1.0 / ( lTwoH + lKpH * lExpmOne );
        lResult[iTime] =
            mMarketZCB.instantaneousForwardRate( mTerms[iTime] ) -
            ( ( lTwoKMean * lExpmOne + lDenom * lFourHHX * lExp ) * lDenom )
                .sum();
    }
    return lResult;
}

}  // namespace MultiFactor
}  // namespace ShortRate
