/**
 * @file interpolate_multi.cpp
 * @brief This implements interpolation function for multi dimensional dataset.
 * @author kakune
 * @date 3/13/2024
 */

#include "math/interpolate_multi.hpp"

namespace Math
{
namespace InterpolateMulti
{

void RBFAbstract::build(
    std::vector<std::shared_ptr<const std::vector<double>>> insRefVars,
    std::shared_ptr<const std::vector<double>> insRefVals )
{
    mNCoeff = insRefVals->size();
    mNDim   = insRefVars.size();
    Math::Vec lValues( mNCoeff );
    for ( std::size_t iPos = 0; iPos < mNCoeff; ++iPos )
    {
        lValues( iPos ) = insRefVals->operator[]( iPos );
    }
    Math::Mat lDistances( mNCoeff, mNCoeff );
    mPoints = std::vector<Math::Vec>( mNCoeff, Math::Vec( mNDim ) );
    for ( std::size_t iPos = 0; iPos < mNCoeff; ++iPos )
    {
        for ( std::size_t iDim = 0; iDim < mNDim; ++iDim )
        {
            mPoints[iPos]( iDim ) = insRefVars[iDim]->operator[]( iPos );
        }
    }
    double lTmpMinDistance = 100.0;
    for ( std::size_t i = 0; i < mNCoeff; ++i )
    {
        for ( std::size_t j = 0; j < i; ++j )
        {
            lDistances( i, j ) = distance( mPoints[i], mPoints[j] );
            lTmpMinDistance = std::min( lTmpMinDistance, lDistances( i, j ) );
        }
    }
    mFactorDistance = mMinDistance / lTmpMinDistance;
    for ( std::size_t i = 0; i < mNCoeff; ++i )
    {
        for ( std::size_t j = 0; j <= i; ++j )
        {
            lDistances( i, j ) = radial( lDistances( i, j ) * mFactorDistance );
            lDistances( j, i ) = lDistances( i, j );
        }
    }
    mCoeffs    = dot( lDistances, lValues );
    lDistances = dot( lDistances, lDistances );
    for ( std::size_t i = 0; i < mNCoeff; ++i )
    {
        lDistances( i, i ) += mFactorDecay;
    }
    solveEqPositiveDefinite( lDistances, mCoeffs );
    mIsBuilt = true;
}
double RBFAbstract::operator()( const std::vector<double>& inVar ) const
{
    Math::Vec lVar( inVar );
    Math::Vec lRadials( mNCoeff );
    for ( std::size_t iPos = 0; iPos < mNCoeff; ++iPos )
    {
        lRadials( iPos ) = radial( distance( lVar, mPoints[iPos] ) );
    }
    return dot( mCoeffs, lRadials );
}
double RBFAbstract::distance( const Math::Vec& inX1,
                              const Math::Vec& inX2 ) const
{
    double lResult = 0.0;
    for ( std::size_t i = 0; i < mNDim; ++i )
    {
        double lTmp = inX1( i ) - inX2( i );
        lResult += lTmp * lTmp;
    }
    return mFactorDistance * lResult;
}
double RBFGaussian::radial( double inDist ) const
{
    return std::exp( -inDist );
}

}  // namespace InterpolateMulti
}  // namespace Math