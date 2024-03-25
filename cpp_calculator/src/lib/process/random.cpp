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

void PathAbstract::initRandomValues()
{
    if ( mRandomValues.size() == mNPath &&
         mRandomValues.at( 0 ).size() == mTerms.size() )
    {
        return;
    }
    mRandomValues.resize( mNPath, std::vector<double>( mTerms.size() ) );
}

void PathAbstract::makePath()
{
    initRandomValues();
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
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
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
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
namespace RandomVec
{

void PathAbstract::initRandomValues()
{
    if ( mRandomValues.size() == mNPath &&
         mRandomValues.at( 0 ).size() == mTerms.size() )
    {
        return;
    }
    mRandomValues.resize(
        mNPath, std::vector<Math::Vec>( mTerms.size(), Math::Vec( mDim ) ) );
}

void PathAbstract::makePath()
{
    initRandomValues();
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
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
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        setIndexTime( iTerm );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mRandomValues[iPath][iTerm] = generateRandomVal();
        }
    }
}

Math::Vec PathBrownPlain::generateRandomVal()
{
    Math::Vec lResult( mDim );
    for ( int i = 0; i < mDim; ++i )
    {
        lResult( i ) = mDistribution( mGenerator ) * mTmpSqrtInterval;
    }
    return lResult;
}

Math::Vec PathBrownAntithetic::generateRandomVal()
{
    if ( !mIsNextNew )
    {
        mIsNextNew = true;
        return -mPrevRandomValue;
    }
    mIsNextNew = false;
    for ( int i = 0; i < mDim; ++i )
    {
        mPrevRandomValue( i ) = mDistribution( mGenerator ) * mTmpSqrtInterval;
    }
    return mPrevRandomValue;
}
}  // namespace RandomVec
}  // namespace Process