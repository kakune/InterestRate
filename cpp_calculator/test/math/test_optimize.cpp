#include <gtest/gtest.h>

#include "math/matrix.hpp"
#include "math/optimize.hpp"

static double testFunctionPoly1( const Math::Vec& inVec )
{
    return 2.0 * inVec[0] * inVec[0] + inVec[1] * inVec[1] +
           inVec[0] * inVec[1] + inVec[0] + inVec[1];
}

double testOptimize( auto inObjectiveFunction, Math::Vec inInitVec,
                     Math::Vec inAnswer )
{
    Math::Vec lDif =
        inAnswer - Math::Optimize::adam( inObjectiveFunction, inInitVec );
    return Math::dot( lDif, lDif );
}

TEST( OptimizeTest, adam )
{
    EXPECT_NEAR( 0.0,
                 testOptimize( testFunctionPoly1, { 0.0, 0.0 },
                               { -1.0 / 7.0, -3.0 / 7.0 } ),
                 1e-6 );
}