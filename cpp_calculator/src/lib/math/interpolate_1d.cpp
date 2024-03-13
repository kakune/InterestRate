/**
 * @file interpolate_1d.cpp
 * @brief This implements interpolation function for 1d.
 * @author kakune
 * @date 1/29/2024
 */

#include "math/interpolate_1d.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

namespace Math
{
namespace Interpolate1D
{

void NewtonSpline::build( std::shared_ptr<const std::vector<double> > insRefXs,
                          std::shared_ptr<const std::vector<double> > insRefYs )
{
    msRefXs = insRefXs;
    msRefYs = insRefYs;
    if ( mNDeg == 0 )
    {
        mIsBuilt = true;
        return;
    }
    std::size_t lSize = insRefXs->size();
    if ( lSize <= mNDeg )
    {
        std::cerr << "Error: Math::Interpolate1D::NewtonSpline::build()"
                  << std::endl
                  << "There are " << lSize
                  << " points, which is not enough to evaluate." << std::endl;
        return;
    }
    std::size_t lMaxIndCoeff = lSize - 1;
    mCoeff.at( 0 ).resize( lSize );
    for ( std::size_t iX = 0; iX < lMaxIndCoeff; ++iX )
    {
        mCoeff.at( 0 ).at( iX ) =
            ( msRefYs->at( iX + 1 ) - msRefYs->at( iX ) ) /
            ( msRefXs->at( iX + 1 ) - msRefXs->at( iX ) );
    }
    for ( std::size_t iDeg = 1; iDeg < mNDeg; ++iDeg )
    {
        mCoeff.at( iDeg ).resize( insRefXs->size() );
        --lMaxIndCoeff;
        for ( std::size_t iX = 0; iX < lMaxIndCoeff; ++iX )
        {
            mCoeff.at( iDeg ).at( iX ) =
                ( mCoeff.at( iDeg - 1 ).at( iX + 1 ) -
                  mCoeff.at( iDeg - 1 ).at( iX ) ) /
                ( msRefXs->at( iX + iDeg + 1 ) - msRefXs->at( iX ) );
        }
        for ( std::size_t iX = lMaxIndCoeff; iX < lSize - 1; ++iX )
        {
            mCoeff.at( iDeg ).at( iX ) = mCoeff.at( iDeg ).at( iX - 1 );
        }
    }
    mIsBuilt = true;
}

double NewtonSpline::operator()( double inX ) const
{
    if ( !mIsBuilt )
    {
        std::cerr << "Error: Math::Interpolate1D::NewtonSpline::operator()"
                  << std::endl
                  << "The spline has NOT been built." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if ( inX < msRefXs->front() || inX > msRefXs->back() )
    {
        std::cerr << "Error: Math::Interpolate1D::NewtonSpline::operator()"
                  << std::endl
                  << "Argument is out of the allowed range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if ( inX == msRefXs->front() ) { return msRefYs->front(); }
    if ( inX == msRefXs->back() ) { return msRefYs->back(); }
    auto lItr = std::upper_bound( msRefXs->begin(), msRefXs->end(), inX );
    --lItr;
    while ( msRefXs->end() - lItr < mNDeg ) { --lItr; }
    std::size_t lInd = lItr - msRefXs->begin();
    double lRes      = msRefYs->at( lInd );
    double lMulDx    = 1.0;
    for ( std::size_t iDx = 0; iDx < mNDeg; ++iDx )
    {
        lMulDx *= inX - *lItr;
        lRes += mCoeff.at( iDx ).at( lInd ) * lMulDx;
        ++lItr;
    }
    return lRes;
}

double NewtonSpline::deriv( double inX, std::size_t inOrder ) const
{
    if ( inOrder == 0 ) { return operator()( inX ); }
    if ( mNDeg == 0 ) { return 0.0; }
    if ( !mIsBuilt )
    {
        std::cerr << "Error: Math::Interpolate1D::NewtonSpline::deriv"
                  << std::endl
                  << "The spline has NOT been built." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if ( inX < msRefXs->front() || inX > msRefXs->back() )
    {
        std::cerr << "Error: Math::Interpolate1D::NewtonSpline::deriv"
                  << std::endl
                  << "Argument is out of the allowed range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto lItr = std::upper_bound( msRefXs->begin(), msRefXs->end(), inX );
    --lItr;
    while ( msRefXs->end() - lItr < mNDeg ) { --lItr; }
    std::size_t lInd = lItr - msRefXs->begin();

    double lRes = 0.0;
    std::vector<double> lMulDxs( inOrder + 1, 0.0 );
    lMulDxs.at( 0 ) = 1.0;
    for ( std::size_t iDx = 0; iDx < mNDeg; ++iDx )
    {
        double lDif = inX - *lItr;
        for ( std::size_t iDeriv = inOrder; iDeriv > 0; --iDeriv )
        {
            lMulDxs.at( iDeriv ) *= lDif;
            lMulDxs.at( iDeriv ) += iDeriv * lMulDxs.at( iDeriv - 1 );
        }
        lMulDxs.at( 0 ) *= lDif;
        lRes += mCoeff.at( iDx ).at( lInd ) * lMulDxs.at( inOrder );
        ++lItr;
    }

    return lRes;
}

void NewtonSpline::buildIntegral()
{
    mCoeffIntegral.resize( mNDeg + 2 );
    for ( std::size_t iDeg = 0; iDeg <= mNDeg + 1; ++iDeg )
    {
        mCoeffIntegral.at( iDeg ).resize( msRefXs->size(), 0.0 );
    }
    for ( std::size_t lIndLeft = 0; lIndLeft < msRefXs->size(); ++lIndLeft )
    {
        if ( msRefXs->size() < mNDeg + lIndLeft ) { break; }
        std::vector<double> lTmpBareCoeff( mNDeg + 1 );
        lTmpBareCoeff.at( 0 )                 = 1.0;
        mCoeffIntegral.at( 0 ).at( lIndLeft ) = msRefYs->at( lIndLeft );
        for ( std::size_t iDeg = 1; iDeg <= mNDeg; ++iDeg )
        {
            double lDif =
                msRefXs->at( lIndLeft ) - msRefXs->at( lIndLeft + iDeg - 1 );
            for ( std::size_t jDeg = iDeg; jDeg > 0; --jDeg )
            {
                lTmpBareCoeff.at( jDeg ) *= lDif;
                lTmpBareCoeff.at( jDeg ) += lTmpBareCoeff.at( jDeg - 1 );
                mCoeffIntegral.at( jDeg ).at( lIndLeft ) +=
                    mCoeff.at( iDeg - 1 ).at( lIndLeft ) *
                    lTmpBareCoeff.at( jDeg );
            }
            lTmpBareCoeff.at( 0 ) = 0;
        }
        for ( std::size_t iDeg = mNDeg + 1; iDeg > 0; --iDeg )
        {
            mCoeffIntegral.at( iDeg ).at( lIndLeft ) =
                mCoeffIntegral.at( iDeg - 1 ).at( lIndLeft ) / double( iDeg );
        }
        mCoeffIntegral.at( 0 ).at( lIndLeft ) = 0;
    }

    mSumIntegral.resize( msRefXs->size(), 0.0 );
    std::size_t lIndLeft = 0;
    for ( std::size_t iX = 1; iX < msRefXs->size(); ++iX )
    {
        mSumIntegral.at( iX ) = mSumIntegral.at( lIndLeft );
        double lDif           = msRefXs->at( iX ) - msRefXs->at( lIndLeft );
        double lPowX          = 1.0;
        for ( std::size_t iDeg = 1; iDeg <= mNDeg + 1; ++iDeg )
        {
            lPowX *= lDif;
            mSumIntegral.at( iX ) +=
                lPowX * mCoeffIntegral.at( iDeg ).at( lIndLeft );
        }
        if ( msRefXs->size() > mNDeg + lIndLeft ) { ++lIndLeft; }
    }
}

double NewtonSpline::integral( double inX ) const
{
    if ( mCoeffIntegral.size() == 0 )
    {
        std::cerr << "Error: Math::Interpolate1D::NewtonSpline::integral"
                  << std::endl
                  << "The coefficients of integral have not calculated."
                  << std::endl
                  << "Run buildIntegral() before integral()." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if ( inX < msRefXs->front() || inX > msRefXs->back() )
    {
        std::cerr << "Error: Math::Interpolate1D::NewtonSpline::integral"
                  << std::endl
                  << "Argument is out of the allowed range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    auto lItr = std::upper_bound( msRefXs->begin(), msRefXs->end(), inX );
    --lItr;
    if ( *lItr == inX ) { return mSumIntegral.at( lItr - msRefXs->begin() ); }
    while ( msRefXs->end() - lItr < mNDeg ) { --lItr; }
    std::size_t lIndLeft = lItr - msRefXs->begin();
    double lResult       = mSumIntegral.at( lIndLeft );
    double lDif          = inX - *lItr;
    double lPowX         = 1.0;
    for ( std::size_t iDeg = 1; iDeg <= mNDeg + 1; ++iDeg )
    {
        lPowX *= lDif;
        lResult += lPowX * mCoeffIntegral.at( iDeg ).at( lIndLeft );
    }
    return lResult;
}

double NewtonSpline::integral( double inA, double inB ) const
{
    return integral( inB ) - integral( inA );
}

}  // namespace Interpolate1D
}  // namespace Math