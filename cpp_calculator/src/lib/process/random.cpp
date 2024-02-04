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

#include "math/interpolate_1d.hpp"

namespace Process
{
namespace Random
{

void PathAbstract::initRandomValues()
{
    if ( mRandomValues.size() == mNPath &&
         mRandomValues.at( 0 ).size() == msTerms->size() )
    {
        return;
    }
    mRandomValues.resize( mNPath, std::vector< double >( msTerms->size() ) );
}

void PathAbstract::makePath()
{
    initRandomValues();
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        setIndexTime( iTerm );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mRandomValues[iPath][iTerm] =
                mRandomValues[iPath][iTerm - 1] + generateRandomVal();
        }
    }
}

void PathAbstract::makeRandomVariables()
{
    initRandomValues();
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        setIndexTime( iTerm );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mRandomValues[iPath][iTerm] = generateRandomVal();
        }
    }
}

double PathBrownPlain::generateRandomVal()
{
    return mDistribution( mGenerator ) * mTmpSqrtInterval;
}

double PathBrownAntithetic::generateRandomVal()
{
    if ( !mIsNextNew )
    {
        mIsNextNew = true;
        return -mPrevRandomValue;
    }
    mIsNextNew       = false;
    mPrevRandomValue = mDistribution( mGenerator ) * mTmpSqrtInterval;
    return mPrevRandomValue;
}

}  // namespace Random
}  // namespace Process