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
#include <memory>
#include <vector>

namespace Math
{
namespace Interpolate1d
{

void NewtonSpline::build( std::shared_ptr<const std::vector<double> > insRefXs,
                          std::shared_ptr<const std::vector<double> > insRefYs )
{
    msRefXs = insRefXs;
    msRefYs = insRefYs;
    if ( mNDeg == 0 ) { return; }
    std::size_t lSize = insRefXs->size();
    if ( lSize <= mNDeg )
    {
        std::cerr << "Error: Math::Interpolate1d::NewtonSpline::build()"
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

double NewtonSpline::operator()( double inX )
{
    if ( !mIsBuilt )
    {
        std::cerr << "Error: Math::Interpolate1d::NewtonSpline::operator()"
                  << std::endl
                  << "The spline has NOT been built." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if ( inX < msRefXs->front() || inX > msRefXs->back() )
    {
        std::cerr << "Error: Math::Interpolate1d::NewtonSpline::operator()"
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

double NewtonSpline::deriv( double inX, std::size_t inOrder )
{
    if ( inOrder == 0 ) { return operator()( inX ); }
    if ( !mIsBuilt )
    {
        std::cerr << "Error: The spline has NOT been built." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if ( inX <= msRefXs->front() || inX >= msRefXs->back() )
    {
        std::cerr << "Error: Argument is out of the allowed range."
                  << std::endl;
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

}  // namespace Interpolate1d
}  // namespace Math