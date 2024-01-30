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

void PathBrownPlain::makePath()
{
    std::random_device lDevice;
    std::mt19937 lGenerator( lDevice() );
    std::normal_distribution< double > lDistribution( 0.0, 1.0 );

    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        double lSqrtInterval =
            std::sqrt( msTerms->at( iTerm ) - msTerms->at( iTerm - 1 ) );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mRandomValues[iPath][iTerm] =
                mRandomValues[iPath][iTerm - 1] +
                lDistribution( lGenerator ) * lSqrtInterval;
        }
    }
}

void PathBrownAntithetic::makePath()
{
    std::random_device lDevice;
    std::mt19937 lGenerator( lDevice() );
    std::normal_distribution< double > lDistribution( 0.0, 1.0 );

    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        double lSqrtInterval =
            std::sqrt( msTerms->at( iTerm ) - msTerms->at( iTerm - 1 ) );
        for ( std::size_t iPath = 0; iPath < mNPath; iPath += 2 )
        {
            mRandomValues[iPath][iTerm] =
                mRandomValues[iPath][iTerm - 1] +
                lDistribution( lGenerator ) * lSqrtInterval;
            if ( iPath < mNPath - 1 )
            {
                mRandomValues[iPath + 1][iTerm] = -mRandomValues[iPath][iTerm];
            }
        }
    }
}

}  // namespace Random
}  // namespace Process