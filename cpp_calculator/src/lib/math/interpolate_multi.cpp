/**
 * @file interpolate_multi.cpp
 * @brief This implements interpolation function for multi dimensional dataset.
 * @author kakune
 * @date 3/13/2024
 */

#include "math/interpolate_multi.hpp"

#include <iostream>
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
    double lTmpMinDistance = 1.0e10;
    for ( std::size_t i = 0; i < mNCoeff; ++i )
    {
        for ( std::size_t j = 0; j < i; ++j )
        {
            lDistances( i, j ) = distance( mPoints[i], mPoints[j] );
            lTmpMinDistance = std::min( lTmpMinDistance, lDistances( i, j ) );
        }
    }
    if ( mMinDistance > 0.0 )
    {
        mFactorDistance = mMinDistance / lTmpMinDistance;
    }
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
double RBFAbstract::deriv( const std::vector<double>& inVar, std::size_t inDim,
                           std::size_t inOrder ) const
{
    if ( inOrder == 0 ) { return operator()( inVar ); }
    if ( inOrder > 2 )
    {
        std::cerr << "Math::InterpolateMulti::RBFAbstract::deriv()" << std::endl
                  << "Derivative of order > 2 is not implemented..."
                  << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    Math::Vec lVar( inVar );
    Math::Vec lRadials( mNCoeff );
    for ( std::size_t iPos = 0; iPos < mNCoeff; ++iPos )
    {
        if ( inOrder == 1 )
        {
            lRadials( iPos ) =
                derivRadial( distance( lVar, mPoints[iPos] ), 1 ) *
                derivDistance( lVar, mPoints[iPos], inDim, 1 );
        }
        if ( inOrder == 2 )
        {
            double lDist     = distance( lVar, mPoints[iPos] );
            double lDDist    = derivDistance( lVar, mPoints[iPos], inDim, 1 );
            lRadials( iPos ) = lDDist * lDDist * derivRadial( lDist, 2 ) +
                               derivDistance( lVar, mPoints[iPos], inDim, 2 ) *
                                   derivRadial( lDist, 1 );
        }
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
double RBFAbstract::derivDistance( const Math::Vec& inX1, const Math::Vec& inX2,
                                   std::size_t inDim,
                                   std::size_t inOrder ) const
{
    if ( inOrder == 0 ) { return distance( inX1, inX2 ); }
    if ( inOrder == 1 )
    {
        return 2.0 * mFactorDistance * ( inX1( inDim ) - inX2( inDim ) );
    }
    if ( inOrder == 2 ) { return 2.0 * mFactorDistance; }
    return 0.0;
}
double RBFGaussian::radial( double inDist ) const
{
    return std::exp( -inDist );
}
double RBFGaussian::derivRadial( double inDist, std::size_t inOrder ) const
{
    if ( inOrder % 2 == 0 ) return radial( inDist );
    return -radial( inDist );
}

}  // namespace InterpolateMulti
}  // namespace Math