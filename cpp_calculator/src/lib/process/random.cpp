/**
 * @file random.cpp
 * @brief This implements random path generators.
 * @author kakune
 * @date 1/29/2024
 */

#include "random.hpp"

#include <cmath>
#include <random>

namespace Process
{
namespace Random
{

double StdBrownPlain::operator()() { return mDistribution( mGenerator ); }

void StdBrownAntithetic::initialize() { mIsNextNew = true; }

double StdBrownAntithetic::operator()()
{
    if ( !mIsNextNew )
    {
        mIsNextNew = true;
        return -mPrevRandomValue;
    }
    mIsNextNew       = false;
    mPrevRandomValue = mDistribution( mGenerator );
    return mPrevRandomValue;
}

}  // namespace Random

namespace RandomVec
{

Math::Vec StdBrownPlain::operator()()
{
    Math::Vec lResult( mDim, 0.0 );
    for ( std::size_t i = mIndStart; i < mDim; ++i )
    {
        lResult( i ) = mDistribution( mGenerator );
    }
    return lResult;
}

void StdBrownAntithetic::initialize( std::size_t inIndStart )
{
    mIndStart  = inIndStart;
    mIsNextNew = true;
}

Math::Vec StdBrownAntithetic::operator()()
{
    if ( !mIsNextNew )
    {
        mIsNextNew = true;
        return -mPrevRandomValue;
    }
    mIsNextNew = false;
    for ( std::size_t i = 0; i < mDim; ++i )
    {
        if ( i < mIndStart ) { mPrevRandomValue( i ) = 0.0; }
        else { mPrevRandomValue( i ) = mDistribution( mGenerator ); }
    }
    return mPrevRandomValue;
}
}  // namespace RandomVec
}  // namespace Process