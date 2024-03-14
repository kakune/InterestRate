#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <memory>

#include "math/interpolate_multi.hpp"

double testRBFGaussianPolynomial( double inX, double inShift, double inDeg,
                                  double inStartX, double inEndX,
                                  std::size_t inNumX )
{
    std::vector<double> lSampleXs( inNumX + 2 );
    double lDifX = ( inEndX - inStartX ) / double( inNumX );
    lSampleXs[0] = inStartX;
    for ( std::size_t iX = 1; iX < inNumX + 2; ++iX )
    {
        lSampleXs[iX] = lSampleXs[iX - 1] + lDifX;
    }
    std::vector<double> lSampleYs( lSampleXs.size() );
    for ( std::size_t iX = 0; iX < lSampleXs.size(); ++iX )
    {
        lSampleYs.at( iX ) = pow( lSampleXs.at( iX ) + inShift, inDeg );
    }
    Math::InterpolateMulti::RBFGaussian lSpline( 0.2, 0.001 );
    lSpline.build( { std::make_shared<std::vector<double>>( lSampleXs ) },
                   std::make_shared<std::vector<double>>( lSampleYs ) );
    std::cout << lSpline( { inX } ) << std::endl;
    return lSpline( { inX } );
}

TEST( InterpolateMultiTest, InterpPolynomial )
{
    EXPECT_NEAR( 27.0,
                 testRBFGaussianPolynomial( 0.0, 3.0, 3.0, -1.0, 1.0, 1001 ),
                 1.0 );
    EXPECT_NEAR( 729.0,
                 testRBFGaussianPolynomial( 6.0, 3.0, 3.0, 0.0, 16.0, 1001 ),
                 1.0 );
    EXPECT_NEAR( 125.0,
                 testRBFGaussianPolynomial( 2.0, 3.0, 3.0, 0.0, 16.0, 1001 ),
                 1.0 );
    EXPECT_NEAR(
        3.0, testRBFGaussianPolynomial( 3.0, 0.0, 1.0, 0.0, 16.0, 1001 ), 1.0 );
    EXPECT_NEAR( 144.0,
                 testRBFGaussianPolynomial( 11.0, 1.0, 2.0, 0.0, 16.0, 1001 ),
                 1.0 );
}