#include <gtest/gtest.h>

#include "math/findroot_1d.hpp"

double testFunctionPoly1( double x ) { return x * x - 2.0; }
double testFunctionPoly2( double x ) { return x * x * x + 5.0; }
double testFunctionPoly3( double x ) { return 0.001 * x + 200.0; }

TEST( FindRoot1DTest, Brent )
{
    EXPECT_NEAR( 1.4142135623730951,
                 Math::FindRoot1D::Brent( testFunctionPoly1, 0.0, 2.0 ), 1e-6 );
    EXPECT_NEAR( -1.709975946676697,
                 Math::FindRoot1D::Brent( testFunctionPoly2, -0.0001, -3000.0 ),
                 1e-6 );
    EXPECT_NEAR(
        -2e5, Math::FindRoot1D::Brent( testFunctionPoly3, -1e10, 1e10 ), 1e-6 );
}