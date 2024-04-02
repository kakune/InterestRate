/**
 * @file random.cpp
 * @brief This implements random path generators.
 * @author kakune
 * @date 1/29/2024
 */

#include "random.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
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
    Math::Vec lResult( mDim );
    for ( int i = 0; i < mDim; ++i )
    {
        lResult( i ) = mDistribution( mGenerator );
    }
    return lResult;
}

void StdBrownAntithetic::initialize() { mIsNextNew = true; }

Math::Vec StdBrownAntithetic::operator()()
{
    if ( !mIsNextNew )
    {
        mIsNextNew = true;
        return -mPrevRandomValue;
    }
    mIsNextNew = false;
    for ( int i = 0; i < mDim; ++i )
    {
        mPrevRandomValue( i ) = mDistribution( mGenerator );
    }
    return mPrevRandomValue;
}
}  // namespace RandomVec
}  // namespace Process