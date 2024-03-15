#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <memory>

#include "math/interpolate_multi.hpp"

double testRBFGaussianPolynomial1D( double inX, double inShift, double inDeg,
                                    double inStartX, double inEndX,
                                    std::size_t inNumX,
                                    std::size_t inDerivOrder = 0 )
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
    double lResult;
    if ( inDerivOrder > 0 )
    {
        lResult = lSpline.deriv( { inX }, 0, inDerivOrder );
    }
    else { lResult = lSpline( { inX } ); }
    std::cout << lResult << std::endl;
    return lResult;
}

double testRBFGaussianPolynomial2D(
    std::pair<double, double> inVar, std::pair<double, double> inShift,
    std::pair<double, double> inDeg, std::pair<double, double> inStart,
    std::pair<double, double> inEnd, std::size_t inNum,
    std::size_t inDerivDim = 0, std::size_t inDerivOrder = 0 )
{
    std::vector<double> lSampleXs, lSampleYs;
    double lDifX = ( inEnd.first - inStart.first ) / double( inNum );
    double lDifY = ( inEnd.second - inStart.second ) / double( inNum );
    double lTmpX = inStart.first;

    for ( std::size_t iX = 0; iX < inNum + 2; ++iX )
    {
        double lTmpY = inStart.second;
        for ( std::size_t iY = 0; iY < inNum + 2; ++iY )
        {
            lSampleXs.push_back( lTmpX );
            lSampleYs.push_back( lTmpY );
            lTmpY += lDifY;
        }
        lTmpX += lDifX;
    }
    std::vector<double> lSampleVals( lSampleXs.size() );
    for ( std::size_t iX = 0; iX < lSampleXs.size(); ++iX )
    {
        lSampleVals.at( iX ) =
            pow( lSampleXs[iX] + inShift.first, inDeg.first ) *
            pow( lSampleYs[iX] + inShift.second, inDeg.second );
    }
    Math::InterpolateMulti::RBFGaussian lSpline( 0.2, 0.001 );
    lSpline.build( { std::make_shared<std::vector<double>>( lSampleXs ),
                     std::make_shared<std::vector<double>>( lSampleYs ) },
                   std::make_shared<std::vector<double>>( lSampleVals ) );
    double lResult;
    if ( inDerivOrder > 0 )
    {
        lResult = lSpline.deriv( { inVar.first, inVar.second }, inDerivDim,
                                 inDerivOrder );
    }
    else { lResult = lSpline( { inVar.first, inVar.second } ); }
    std::cout << lResult << std::endl;
    return lResult;
}

TEST( InterpolateMultiTest, InterpPolynomial1D )
{
    // (x+3)^3 at x=0
    EXPECT_NEAR( 27.0,
                 testRBFGaussianPolynomial1D( 0.0, 3.0, 3.0, -1.0, 1.0, 100 ),
                 1.0 );
    EXPECT_NEAR(
        27.0, testRBFGaussianPolynomial1D( 0.0, 3.0, 3.0, -1.0, 1.0, 100, 1 ),
        1.0 );
    EXPECT_NEAR(
        18.0, testRBFGaussianPolynomial1D( 0.0, 3.0, 3.0, -1.0, 1.0, 100, 2 ),
        1.0 );

    // (x+3)^3 at x=6
    EXPECT_NEAR( 729.0,
                 testRBFGaussianPolynomial1D( 6.0, 3.0, 3.0, 0.0, 16.0, 100 ),
                 1.0 );
    EXPECT_NEAR(
        243.0, testRBFGaussianPolynomial1D( 6.0, 3.0, 3.0, 0.0, 16.0, 100, 1 ),
        1.0 );
    EXPECT_NEAR(
        54.0, testRBFGaussianPolynomial1D( 6.0, 3.0, 3.0, 0.0, 16.0, 100, 2 ),
        1.0 );

    // (x+1)^8 at x=1
    EXPECT_NEAR( 256.0,
                 testRBFGaussianPolynomial1D( 1.0, 1.0, 8.0, 0.0, 16.0, 100 ),
                 1.0 );
    EXPECT_NEAR(
        1024.0, testRBFGaussianPolynomial1D( 1.0, 1.0, 8.0, 0.0, 16.0, 100, 1 ),
        1.0 );
    EXPECT_NEAR(
        3584.0, testRBFGaussianPolynomial1D( 1.0, 1.0, 8.0, 0.0, 16.0, 100, 2 ),
        1.0 );

    // x at x=3
    EXPECT_NEAR( 3.0,
                 testRBFGaussianPolynomial1D( 3.0, 0.0, 1.0, 0.0, 16.0, 100 ),
                 1.0 );
    EXPECT_NEAR(
        1.0, testRBFGaussianPolynomial1D( 3.0, 0.0, 1.0, 0.0, 16.0, 100, 1 ),
        1.0 );
    EXPECT_NEAR(
        0.0, testRBFGaussianPolynomial1D( 3.0, 0.0, 1.0, 0.0, 16.0, 100, 2 ),
        1.0 );

    // (x+1)^2 at x=11
    EXPECT_NEAR(
        144.0, testRBFGaussianPolynomial1D( 11.0, 1.0, 2.0, 0.0, 1000.0, 500 ),
        1.0 );
    EXPECT_NEAR(
        24.0, testRBFGaussianPolynomial1D( 11.0, 1.0, 2.0, 0.0, 100.0, 500, 1 ),
        1.0 );
    EXPECT_NEAR(
        2.0, testRBFGaussianPolynomial1D( 11.0, 1.0, 2.0, 0.0, 100.0, 500, 2 ),
        1.0 );

    // 1 at x=11
    EXPECT_NEAR(
        1.0, testRBFGaussianPolynomial1D( 11.0, 1.0, 0.0, 0.0, 1000.0, 100 ),
        1.0 );
    EXPECT_NEAR(
        0.0, testRBFGaussianPolynomial1D( 11.0, 1.0, 0.0, 0.0, 100.0, 100, 1 ),
        1.0 );
    EXPECT_NEAR(
        0.0, testRBFGaussianPolynomial1D( 11.0, 1.0, 0.0, 0.0, 100.0, 100, 2 ),
        1.0 );
}

TEST( InterpolateMultiTest, InterpPolynomial2D )
{
    //(x+3)^3 (y+3)^3 at (0,0)
    EXPECT_NEAR(
        729.0,
        testRBFGaussianPolynomial2D( { 0.0, 0.0 }, { 3.0, 3.0 }, { 3.0, 3.0 },
                                     { -1.0, -1.0 }, { 1.0, 1.0 }, 30 ),
        1.0 );
    EXPECT_NEAR(
        729.0,
        testRBFGaussianPolynomial2D( { 0.0, 0.0 }, { 3.0, 3.0 }, { 3.0, 3.0 },
                                     { -1.0, -1.0 }, { 1.0, 1.0 }, 50, 0, 1 ),
        1.0 );
    EXPECT_NEAR(
        486.0,
        testRBFGaussianPolynomial2D( { 0.0, 0.0 }, { 3.0, 3.0 }, { 3.0, 3.0 },
                                     { -1.0, -1.0 }, { 1.0, 1.0 }, 50, 1, 2 ),
        24.0 );

    //(x+1) at (0,0)
    EXPECT_NEAR(
        1.0,
        testRBFGaussianPolynomial2D( { 0.0, 0.0 }, { 1.0, 2.0 }, { 1.0, 0.0 },
                                     { -1.0, -1.0 }, { 1.0, 1.0 }, 30 ),
        0.01 );
    EXPECT_NEAR(
        1.0,
        testRBFGaussianPolynomial2D( { 0.0, 0.0 }, { 1.0, 2.0 }, { 1.0, 0.0 },
                                     { -1.0, -1.0 }, { 1.0, 1.0 }, 50, 0, 1 ),
        0.02 );
    EXPECT_NEAR(
        0.0,
        testRBFGaussianPolynomial2D( { 0.0, 0.0 }, { 1.0, 2.0 }, { 1.0, 0.0 },
                                     { -1.0, -1.0 }, { 1.0, 1.0 }, 50, 1, 2 ),
        0.05 );

    // x y^2 at (2,3)
    EXPECT_NEAR(
        18.0,
        testRBFGaussianPolynomial2D( { 2.0, 3.0 }, { 0.0, 0.0 }, { 1.0, 2.0 },
                                     { -1.0, -1.0 }, { 10.0, 10.0 }, 30 ),
        0.18 );
    EXPECT_NEAR(
        9.0,
        testRBFGaussianPolynomial2D( { 2.0, 3.0 }, { 0.0, 0.0 }, { 1.0, 2.0 },
                                     { -1.0, -1.0 }, { 10.0, 10.0 }, 50, 0, 1 ),
        0.36 );
    EXPECT_NEAR(
        4.0,
        testRBFGaussianPolynomial2D( { 2.0, 3.0 }, { 0.0, 0.0 }, { 1.0, 2.0 },
                                     { -1.0, -1.0 }, { 10.0, 10.0 }, 50, 1, 2 ),
        0.2 );

    // 1 at (4,1)
    EXPECT_NEAR(
        1.0,
        testRBFGaussianPolynomial2D( { 4.0, 1.0 }, { 0.0, 0.0 }, { 0.0, 0.0 },
                                     { -1.0, -1.0 }, { 10.0, 10.0 }, 30 ),
        0.1 );
    EXPECT_NEAR(
        0.0,
        testRBFGaussianPolynomial2D( { 4.0, 1.0 }, { 0.0, 0.0 }, { 0.0, 0.0 },
                                     { -1.0, -1.0 }, { 10.0, 10.0 }, 30, 0, 1 ),
        0.1 );
    EXPECT_NEAR(
        0.0,
        testRBFGaussianPolynomial2D( { 4.0, 1.0 }, { 0.0, 0.0 }, { 0.0, 0.0 },
                                     { -1.0, -1.0 }, { 10.0, 10.0 }, 30, 1, 2 ),
        0.1 );
}