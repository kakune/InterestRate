#include <gtest/gtest.h>

#include <cmath>

#include "math/integral_1d.hpp"

static double batteryTest1( double x ) { return exp( x ); }
static double batteryTest2( double x ) { return x > 0.3 ? 1 : 0; }
static double batteryTest3( double x ) { return sqrt( x ); }
static double batteryTest4( double x )
{
    return 23.0 * cosh( x ) / 25.0 - cos( x );
}
static double batteryTest5( double x )
{
    return 1.0 / ( x * x * x * x + x * x + 0.9 );
}
static double batteryTest6( double x ) { return x * sqrt( x ); }
static double batteryTest7( double x ) { return 1.0 / sqrt( x ); }
static double batteryTest8( double x ) { return 1.0 / ( 1.0 + x * x * x * x ); }
static double batteryTest9( double x )
{
    return 2.0 / ( 2.0 + sin( 10 * M_PI * x ) );
}
static double batteryTest10( double x ) { return 1.0 / ( 1.0 + x ); }
static double batteryTest11( double x ) { return 1.0 / ( 1.0 + exp( x ) ); }
static double batteryTest12( double x )
{
    return ( x == 0.0 ) ? 1 : x / ( exp( x ) - 1 );
}
static double batteryTest13( double x )
{
    return ( x == 0.0 ) ? 100.0 : sin( 100.0 * M_PI * x ) / ( M_PI * x );
}
static double batteryTest14( double x )
{
    return sqrt( 50.0 ) * exp( -50.0 * M_PI * x * x );
}
static double batteryTest15( double x ) { return 25.0 * exp( -25.0 * x ); }
static double batteryTest16( double x )
{
    return 50.0 / ( M_PI * ( 2500.0 * x * x + 1 ) );
}
static double batteryTest17( double x )
{
    if ( x == 0.0 ) return 50.0;
    return 50.0 * ( sin( 50.0 * M_PI * x ) / ( 50.0 * M_PI * x ) ) *
           ( sin( 50.0 * M_PI * x ) / ( 50.0 * M_PI * x ) );
}
static double batteryTest18( double x )
{
    return cos( cos( x ) + 3.0 * sin( x ) + 2.0 * cos( 2.0 * x ) +
                3.0 * sin( 2.0 * x ) + 3.0 * cos( 3.0 * x ) );
}
static double batteryTest19( double x ) { return ( x > 1e-15 ) ? log( x ) : 0; }
static double batteryTest20( double x ) { return 1.0 / ( 1.005 + x * x ); }
static double batteryTest21( double x )
{
    double lRes = 0.0;
    for ( int i = 1; i <= 3; ++i )
    {
        lRes += 1.0 / cosh( pow( 20.0, i ) * ( x - 0.2 * i ) );
    }
    return lRes;
}
static double batteryTest22( double x )
{
    return 4.0 * M_PI * M_PI * x * sin( 20.0 * M_PI * x ) *
           cos( 2.0 * M_PI * x );
}
static double batteryTest23( double x )
{
    return 1.0 / ( 1.0 + ( 230.0 * x - 30.0 ) * ( 230.0 * x - 30.0 ) );
}

static double lRes1 = exp( 1.0 ) - 1.0;
static double lRes2 = 0.7;
static double lRes3 = 2.0 / 3.0;
static double lRes4 = -2.0 * sin( 1 ) + 46.0 * sinh( 1 ) / 25.0;
static double lRes5 = 1.582232963729673;
static double lRes6 = 2.0 / 5.0;
static double lRes7 = 2.0;
static double lRes8 =
    ( M_PI + 2.0 * acosh( sqrt( 2.0 ) ) ) / ( 4.0 * sqrt( 2.0 ) );
static double lRes9  = 2.0 / sqrt( 3.0 );
static double lRes10 = log( 2.0 );
static double lRes11 = 1.0 + log( 2.0 ) - log( 1.0 + M_E );
static double lRes12 = 0.7775046341122478;
static double lRes13 = 0.00909863753916684;
static double lRes14 = 0.5;
static double lRes15 = 1.0;
static double lRes16 = atan( 500 ) / M_PI;
static double lRes17 = 0.1121393037416375;
static double lRes18 = 0.838676342694364;
static double lRes19 = -1.0;
static double lRes20 =
    20.0 * sqrt( 2.0 / 201.0 ) * atan( 10.0 * sqrt( 2.0 / 201.0 ) );
static double lRes21 = 0.1634949430186372;
static double lRes22 = -20.0 * M_PI / 99.0;
static double lRes23 = ( atan( 200.0 ) + atan( 30.0 ) ) / 230.0;

static double lIRes5 =
    ( M_PI / 3.0 ) * sqrt( 10.0 * ( 3.0 * sqrt( 10.0 ) - 5.0 ) / 13.0 );
static double lIRes8  = M_PI / sqrt( 2.0 );
static double lIRes13 = 1.0;
static double lIRes14 = 1.0;
static double lIRes16 = 1.0;
static double lIRes17 = 1.0;
static double lIRes20 = 10.0 * M_PI * sqrt( 2.0 / 201.0 );
static double lIRes23 = M_PI / 230.0;

static double lURes5  = 1.763856757494504;
static double lURes8  = M_PI / sqrt( 8.0 );
static double lURes11 = log( 2.0 );
static double lURes12 = M_PI * M_PI / 6.0;
static double lURes13 = 0.0101118;
static double lURes14 = 0.5;
static double lURes15 = 1.0;
static double lURes16 = 0.5;
static double lURes17 = 0.11315249504859;
static double lURes20 = 5.0 * sqrt( 2.0 / 201.0 ) *
                        ( M_PI + 2.0 * atan( 10.0 * sqrt( 2.0 / 201.0 ) ) );
static double lURes23 = ( M_PI + 2.0 * atan( 30.0 ) ) / 460.0;

static double gTolAbs = 1e-12;
static double gTolRel = 1e-14;
static double gEps    = std::numeric_limits<double>::epsilon();
static double dblAdaptFinite( auto inFunc, double inMin, double inMax )
{
    return Math::Integral::FiniteInterval::doublyAdaptiveNewtonCotes(
        inFunc, inMin, inMax, gTolAbs, gTolRel );
}
static double deFinite( auto inFunc, double inMin, double inMax )
{
    return Math::Integral::FiniteInterval::DEFormula( inFunc, inMin, inMax,
                                                      gTolRel );
}

TEST( Integral1DTest, DE )
{
    EXPECT_NEAR( lRes1, deFinite( batteryTest1, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes2, deFinite( batteryTest2, 0.0, 1.0 ), 1e-7 );
    EXPECT_NEAR( lRes3, deFinite( batteryTest3, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes4, deFinite( batteryTest4, -1.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes5, deFinite( batteryTest5, -1.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes6, deFinite( batteryTest6, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes7, deFinite( batteryTest7, gEps, 1.0 ), 1e-7 );
    EXPECT_NEAR( lRes8, deFinite( batteryTest8, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes9, deFinite( batteryTest9, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes10, deFinite( batteryTest10, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes11, deFinite( batteryTest11, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes12, deFinite( batteryTest12, gEps, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes13, deFinite( batteryTest13, 0.1, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes14, deFinite( batteryTest14, 0.0, 10.0 ), gTolAbs );
    EXPECT_NEAR( lRes15, deFinite( batteryTest15, 0.0, 10.0 ), gTolAbs );
    EXPECT_NEAR( lRes16, deFinite( batteryTest16, 0.0, 10.0 ), gTolAbs );
    EXPECT_NEAR( lRes17, deFinite( batteryTest17, 0.01, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes18, deFinite( batteryTest18, 0.0, M_PI ), gTolAbs );
    EXPECT_NEAR( lRes19, deFinite( batteryTest19, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes20, deFinite( batteryTest20, -1.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes21, deFinite( batteryTest21, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes22, deFinite( batteryTest22, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes23, deFinite( batteryTest23, 0.0, 1.0 ), gTolAbs );
}

TEST( Integral1DTest, DoublyAdaptive )
{
    EXPECT_NEAR( lRes1, dblAdaptFinite( batteryTest1, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes2, dblAdaptFinite( batteryTest2, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes3, dblAdaptFinite( batteryTest3, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes4, dblAdaptFinite( batteryTest4, -1.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes5, dblAdaptFinite( batteryTest5, -1.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes6, dblAdaptFinite( batteryTest6, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes7, dblAdaptFinite( batteryTest7, gEps, 1.0 ), 1e-8 );
    EXPECT_NEAR( lRes8, dblAdaptFinite( batteryTest8, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes9, dblAdaptFinite( batteryTest9, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes10, dblAdaptFinite( batteryTest10, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes11, dblAdaptFinite( batteryTest11, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes12, dblAdaptFinite( batteryTest12, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes13, dblAdaptFinite( batteryTest13, 0.1, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes14, dblAdaptFinite( batteryTest14, 0.0, 10.0 ), gTolAbs );
    EXPECT_NEAR( lRes15, dblAdaptFinite( batteryTest15, 0.0, 10.0 ), gTolAbs );
    EXPECT_NEAR( lRes16, dblAdaptFinite( batteryTest16, 0.0, 10.0 ), gTolAbs );
    EXPECT_NEAR( lRes17, dblAdaptFinite( batteryTest17, 0.01, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes18, dblAdaptFinite( batteryTest18, 0.0, M_PI ), gTolAbs );
    EXPECT_NEAR( lRes19, dblAdaptFinite( batteryTest19, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes20, dblAdaptFinite( batteryTest20, -1.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes21, dblAdaptFinite( batteryTest21, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes22, dblAdaptFinite( batteryTest22, 0.0, 1.0 ), gTolAbs );
    EXPECT_NEAR( lRes23, dblAdaptFinite( batteryTest23, 0.0, 1.0 ), gTolAbs );
}

static double dblAdaptInfinite( auto inFunc )
{
    return Math::Integral::InfiniteInterval::doublyAdaptiveNewtonCotes(
        inFunc, gTolAbs, gTolRel );
}
static double deInfinite( auto inFunc )
{
    return Math::Integral::InfiniteInterval::DEFormula( inFunc, gTolRel );
}

TEST( Integral1DTest, DEInfinite )
{
    EXPECT_NEAR( lIRes5, deInfinite( batteryTest5 ), gTolAbs );
    EXPECT_NEAR( lIRes8, deInfinite( batteryTest8 ), gTolAbs );
    EXPECT_NEAR( lIRes13, deInfinite( batteryTest13 ), 1e-2 );
    EXPECT_NEAR( lIRes14, deInfinite( batteryTest14 ), gTolAbs );
    EXPECT_NEAR( lIRes16, deInfinite( batteryTest16 ), gTolAbs );
    EXPECT_NEAR( lIRes17, deInfinite( batteryTest17 ), 1e-10 );
    EXPECT_NEAR( lIRes20, deInfinite( batteryTest20 ), gTolAbs );
    EXPECT_NEAR( lIRes23, deInfinite( batteryTest23 ), gTolAbs );
}

TEST( Integral1DTest, DoublyAdaptiveInfinite )
{
    EXPECT_NEAR( lIRes5, dblAdaptInfinite( batteryTest5 ), gTolAbs );
    EXPECT_NEAR( lIRes8, dblAdaptInfinite( batteryTest8 ), gTolAbs );
    // EXPECT_NEAR( lIRes13, dblAdaptInfinite( batteryTest13 ), gTolAbs );
    EXPECT_NEAR( lIRes14, dblAdaptInfinite( batteryTest14 ), gTolAbs );
    EXPECT_NEAR( lIRes16, dblAdaptInfinite( batteryTest16 ), gTolAbs );
    EXPECT_NEAR( lIRes17, dblAdaptInfinite( batteryTest17 ), 1e-4 );
    EXPECT_NEAR( lIRes20, dblAdaptInfinite( batteryTest20 ), gTolAbs );
    EXPECT_NEAR( lIRes23, dblAdaptInfinite( batteryTest23 ), gTolAbs );
}

static double dblAdaptUInfinite( auto inFunc, double inMin )
{
    return Math::Integral::UpperInfiniteInterval::doublyAdaptiveNewtonCotes(
        inFunc, inMin, gTolAbs, gTolRel );
}
static double deUInfinite( auto inFunc, double inMin )
{
    return Math::Integral::UpperInfiniteInterval::DEFormula( inFunc, inMin,
                                                             gTolRel );
}

TEST( Integral1DTest, DEUInfinite )
{
    EXPECT_NEAR( lURes5, deUInfinite( batteryTest5, -1.0 ), gTolAbs );
    EXPECT_NEAR( lURes8, deUInfinite( batteryTest8, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes11, deUInfinite( batteryTest11, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes12, deUInfinite( batteryTest12, gEps ), 1e-7 );
    // EXPECT_NEAR( lURes13, deUInfinite( batteryTest13, 0.1 ), gTolAbs );
    EXPECT_NEAR( lURes14, deUInfinite( batteryTest14, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes15, deUInfinite( batteryTest15, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes16, deUInfinite( batteryTest16, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes17, deUInfinite( batteryTest17, 0.01 ), 1e-10 );
    EXPECT_NEAR( lURes20, deUInfinite( batteryTest20, -1.0 ), gTolAbs );
    EXPECT_NEAR( lURes23, deUInfinite( batteryTest23, 0.0 ), gTolAbs );
}
TEST( Integral1DTest, DoublyAdaptiveUInfinite )
{
    EXPECT_NEAR( lURes5, dblAdaptUInfinite( batteryTest5, -1.0 ), gTolAbs );
    EXPECT_NEAR( lURes8, dblAdaptUInfinite( batteryTest8, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes11, dblAdaptUInfinite( batteryTest11, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes12, dblAdaptUInfinite( batteryTest12, 0.0 ), gTolAbs );
    // EXPECT_NEAR( lURes13, dblAdaptUInfinite( batteryTest13, 0.1 ), gTolAbs );
    EXPECT_NEAR( lURes14, dblAdaptUInfinite( batteryTest14, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes15, dblAdaptUInfinite( batteryTest15, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes16, dblAdaptUInfinite( batteryTest16, 0.0 ), gTolAbs );
    EXPECT_NEAR( lURes17, dblAdaptUInfinite( batteryTest17, 0.01 ), 1e-4 );
    EXPECT_NEAR( lURes20, dblAdaptUInfinite( batteryTest20, -1.0 ), gTolAbs );
    EXPECT_NEAR( lURes23, dblAdaptUInfinite( batteryTest23, 0.0 ), gTolAbs );
}