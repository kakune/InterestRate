#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <memory>

#include "math/interpolate_1d.hpp"

double testInterpPolynomial( double inX, double inShift, double inDeg,
                             std::size_t inInterpDeg,
                             std::vector<double> inSampleXs )
{
    std::vector<double> lSampleYs( inSampleXs.size() );
    for ( std::size_t iX = 0; iX < inSampleXs.size(); ++iX )
    {
        lSampleYs.at( iX ) = pow( inSampleXs.at( iX ) + inShift, inDeg );
    }
    Math::Interpolate1D::NewtonSpline lSpline( inInterpDeg );
    lSpline.build( std::make_shared<std::vector<double>>( inSampleXs ),
                   std::make_shared<std::vector<double>>( lSampleYs ) );
    return lSpline( inX );
}

double testDerivPolynomial( double inX, double inShift, double inDeg,
                            std::size_t inInterpDeg, std::size_t inOrderDeriv,
                            std::vector<double> inSampleXs )
{
    std::vector<double> lSampleYs( inSampleXs.size() );
    for ( std::size_t iX = 0; iX < inSampleXs.size(); ++iX )
    {
        lSampleYs.at( iX ) = pow( inSampleXs.at( iX ) + inShift, inDeg );
    }
    Math::Interpolate1D::NewtonSpline lSpline( inInterpDeg );
    lSpline.build( std::make_shared<std::vector<double>>( inSampleXs ),
                   std::make_shared<std::vector<double>>( lSampleYs ) );
    return lSpline.deriv( inX, inOrderDeriv );
}

double testIntegralPolynomial( double inLeftInterval, double inRightInterval,
                               double inShift, double inDeg,
                               std::size_t inInterpDeg,
                               std::vector<double> inSampleXs )
{
    std::vector<double> lSampleYs( inSampleXs.size() );
    for ( std::size_t iX = 0; iX < inSampleXs.size(); ++iX )
    {
        lSampleYs.at( iX ) = pow( inSampleXs.at( iX ) + inShift, inDeg );
    }
    Math::Interpolate1D::NewtonSpline lSpline( inInterpDeg );
    lSpline.build( std::make_shared<std::vector<double>>( inSampleXs ),
                   std::make_shared<std::vector<double>>( lSampleYs ) );
    lSpline.buildIntegral();
    return lSpline.integral( inLeftInterval, inRightInterval );
}

TEST( Interpolate1DTest, InterpPolynomial )
{
    EXPECT_NEAR(
        125.0, testInterpPolynomial( 2.0, 3.0, 3.0, 3, { 0.0, 2.0, 4.0, 6.0 } ),
        1e-6 );
    EXPECT_NEAR(
        2.0, testInterpPolynomial( 3.0, 0.0, 1.0, 0, { 0.0, 2.0, 4.0, 6.0 } ),
        1e-6 );
    EXPECT_NEAR( 144.0,
                 testInterpPolynomial(
                     11.0, 1.0, 2.0, 6,
                     { 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0 } ),
                 1e-6 );
}

TEST( Interpolate1DTest, DerivPolynomial )
{
    EXPECT_NEAR(
        75.0,
        testDerivPolynomial( 2.0, 3.0, 3.0, 3, 1, { 0.0, 2.0, 4.0, 6.0 } ),
        1e-6 );
    EXPECT_NEAR(
        0.0, testDerivPolynomial( 3.0, 0.0, 1.0, 0, 1, { 0.0, 2.0, 4.0, 6.0 } ),
        1e-6 );
    EXPECT_NEAR( 2.0,
                 testDerivPolynomial(
                     11.0, 1.0, 2.0, 6, 2,
                     { 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0 } ),
                 1e-6 );
}

TEST( Interpolate1DTest, IntegralPolynomial )
{
    EXPECT_NEAR(
        536.25,
        testIntegralPolynomial( 1.0, 4.0, 3.0, 3.0, 3, { 0.0, 2.0, 4.0, 6.0 } ),
        1e-6 );
    EXPECT_NEAR(
        12.0,
        testIntegralPolynomial( 0.0, 6.0, 0.0, 1.0, 0, { 0.0, 2.0, 4.0, 6.0 } ),
        1e-6 );
    EXPECT_NEAR( 242.6666666666667,
                 testIntegralPolynomial(
                     9.0, 11.0, 1.0, 2.0, 6,
                     { 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0 } ),
                 1e-6 );
}