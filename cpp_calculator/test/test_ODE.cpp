#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "math/ODE.hpp"

// the result is $C e^{e^{x}-x}$
double firstOrderODEForTest( double inX, double inY )
{
    return inY * ( std::exp( inX ) - 1 );
}

// the result is ${C cos(x^2/2) + D sin(x^2/2), D cos(x^2/2) - C sin(x^2/2)}$
std::vector<double> SIMLFirstOrderODEForTest( double inX,
                                              std::vector<double> inY )
{
    std::vector<double> lResult( 2 );
    lResult.at( 0 ) = inX * inY.at( 1 );
    lResult.at( 1 ) = -inX * inY.at( 0 );
    return lResult;
}

// the result is $ C e^x + D xe^x + cos(x)/2 $
double secondOrderODEForTest( double inX, double inY, double inDY )
{
    return 2 * inDY - inY + std::sin( inX );
}

// the result is ${2e^x + A + x(B+Cx) - cos(x), e^x + Cx + D}$
std::vector<double> SIMLSecondOrderODEForTest( double inX,
                                               std::vector<double> inY,
                                               std::vector<double> inDY )
{
    std::vector<double> lResults( 2 );
    lResults.at( 0 ) = 2.0 * inDY.at( 1 ) + std::cos( inX );
    lResults.at( 1 ) = std::exp( inX );
    return lResults;
}

double testFirstOrderODE( double inInitX, double inInitY, double inEndX,
                          double inRelTol, double inMaxDif = 1.0 )
{
    auto [lResXs, lResYs] = Math::ODE::solveRungeKutta45<firstOrderODEForTest>(
        inInitX, inInitY, inEndX, inRelTol, inMaxDif );
    // for ( std::size_t i = 0; i < lResXs.size(); ++i )
    // {
    //     std::cout << lResXs.at( i ) << "," << lResYs.at( i ) << std::endl;
    // }
    return lResYs.back();
}

std::vector<double> testSIMLFirstOrderODE( double inInitX,
                                           std::vector<double> inInitY,
                                           double inEndX, double inRelTol,
                                           double inMaxDif = 1.0 )
{
    auto lResults = Math::ODE::solveSIMLRungeKutta45<SIMLFirstOrderODEForTest>(
        inInitX, inInitY, inEndX, inRelTol, inMaxDif );
    std::vector<double> lAnswer;
    // for ( std::size_t i = 0; i < lResults.at( 0 ).size(); ++i )
    // {
    //     for ( std::size_t j = 0; j < lResults.size(); ++j )
    //     {
    //         if ( j > 0 ) { std::cout << ","; }
    //         std::cout << lResults.at( j ).at( i );
    //     }
    //     std::cout << std::endl;
    // }
    for ( std::size_t j = 1; j < lResults.size(); ++j )
    {
        lAnswer.push_back( lResults.at( j ).back() );
    }

    return lAnswer;
}

double testSecondOrderODE( double inInitX, double inInitY, double inInitDY,
                           double inEndX, double inRelTol,
                           double inMaxDif = 1.0 )
{
    auto [lResXs, lResYs, lResDYs] =
        Math::ODE::solveSecondOrderRungeKutta45<secondOrderODEForTest>(
            inInitX, inInitY, inInitDY, inEndX, inRelTol, inMaxDif );
    // for ( std::size_t i = 0; i < lResXs.size(); ++i )
    // {
    //     std::cout << lResXs.at( i ) << "," << lResYs.at( i ) << std::endl;
    // }
    return lResYs.back();
}

std::vector<double> testSIMLSecondOrderODE( double inInitX,
                                            std::vector<double> inInitY,
                                            std::vector<double> inInitDY,
                                            double inEndX, double inRelTol,
                                            double inMaxDif = 1.0 )
{
    auto lResults =
        Math::ODE::solveSIMLSecondOrderRungeKutta45<SIMLSecondOrderODEForTest>(
            inInitX, inInitY, inInitDY, inEndX, inRelTol, inMaxDif );
    std::vector<double> lAnswer;
    // for ( std::size_t i = 0; i < lResults.at( 0 ).size(); ++i )
    // {
    //     for ( std::size_t j = 0; j < lResults.size(); ++j )
    //     {
    //         if ( j > 0 ) { std::cout << ","; }
    //         std::cout << lResults.at( j ).at( i );
    //     }
    //     std::cout << std::endl;
    // }
    for ( std::size_t j = 1; j < lResults.size(); ++j )
    {
        lAnswer.push_back( lResults.at( j ).back() );
    }

    return lAnswer;
}

TEST( ODETest, FirstOrder )
{
    EXPECT_NEAR( 2.050906372692501, testFirstOrderODE( 0.0, 1.0, 1.0, 1e-6 ),
                 1e-5 );
    EXPECT_NEAR( 1.444667861009766, testFirstOrderODE( 0.0, 1.0, -1.0, 1e-6 ),
                 1e-5 );
    EXPECT_NEAR( 579.1374253149154,
                 testFirstOrderODE( 10.0, 64.0, 10.0001, 1e-6 ), 1e-3 );
    EXPECT_NEAR( 1.759637101140221e-8,
                 testFirstOrderODE( 10.0, 64.0, 9.999, 1e-6 ), 1e-10 );
}

TEST( ODETest, SIMLFirstOrder )
{
    EXPECT_NEAR( 0.479425538604203,
                 testSIMLFirstOrderODE( 0.0, { 0.0, 1.0 }, 1.0, 1e-6 ).at( 0 ),
                 1e-6 );
    EXPECT_NEAR( 0.8775825618903728,
                 testSIMLFirstOrderODE( 0.0, { 0.0, 1.0 }, 1.0, 1e-6 ).at( 1 ),
                 1e-6 );
    EXPECT_NEAR(
        0.4926006891650204,
        testSIMLFirstOrderODE( -10.0, { 0.5, -0.2 }, -20.0, 1e-6 ).at( 0 ),
        1e-6 );
    EXPECT_NEAR(
        0.2175880535189073,
        testSIMLFirstOrderODE( -10.0, { 0.5, -0.2 }, -20.0, 1e-6 ).at( 1 ),
        1e-6 );
}

TEST( ODETest, SecondOrder )
{
    EXPECT_NEAR( 0.2701511529340699,
                 testSecondOrderODE( 0.0, 0.0, 0.0, 1.0, 1e-6 ), 1e-5 );
    EXPECT_NEAR( -0.0977282882373725,
                 testSecondOrderODE( 0.0, 0.0, 0.0, -1.0, 1e-6 ), 1e-5 );
    EXPECT_NEAR( 0.1761512469178041,
                 testSecondOrderODE( 2.5, 1.0, -5.0, -7.5, 1e-6 ), 1e-5 );
    EXPECT_NEAR( -4319.555952218765,
                 testSecondOrderODE( 2.5, 1.0, -5.0, 7.5, 1e-6 ), 1e-3 );
}

TEST( ODETest, SIMLSecondOrder )
{
    EXPECT_NEAR(
        0.896261351049951,
        testSIMLSecondOrderODE( 0.0, { 0.0, 0.0 }, { 0.0, 0.0 }, 1.0, 1e-6 )
            .at( 0 ),
        1e-6 );
    EXPECT_NEAR(
        0.7182818284590451,
        testSIMLSecondOrderODE( 0.0, { 0.0, 0.0 }, { 0.0, 0.0 }, 1.0, 1e-6 )
            .at( 1 ),
        1e-6 );
    EXPECT_NEAR(
        -456.3855459491489,
        testSIMLSecondOrderODE( 2.0, { 6.0, -4.0 }, { -2.0, -8.0 }, -4.0, 1e-6 )
            .at( 0 ),
        1e-3 );
    EXPECT_NEAR(
        80.963596133542,
        testSIMLSecondOrderODE( 2.0, { 6.0, -4.0 }, { -2.0, -8.0 }, -4.0, 1e-6 )
            .at( 1 ),
        1e-3 );
}
