/**
 * @file ODE.tpp
 * @brief This implements 1-dim ode solver.
 * @author kakune
 * @date 3/11/2024
 */

#ifndef MATH_ODE_1D_TPP
#define MATH_ODE_1D_TPP

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "math/ODE.hpp"

namespace Math
{
namespace ODE
{

static constexpr double gCoeffKX[6]    = { 0.0,         0.25, 3.0 / 8.0,
                                           12.0 / 13.0, 1.0,  0.5 };
static constexpr double gCoeffKY[6][5] = {
    { 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 0.25, 0.0, 0.0, 0.0, 0.0 },
    { 3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0 },
    { 1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0 },
    { 439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0 },
    { -8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0 } };

static constexpr double gCoeffY4[6] = { 25.0 / 216.0,    0.0,  1408.0 / 2565.0,
                                        2197.0 / 4104.0, -0.2, 0.0 };
static constexpr double gCoeffY5[6] = { 16.0 / 135.0,     0.0,
                                        6656.0 / 12825.0, 28561.0 / 56430.0,
                                        -9.0 / 50.0,      2.0 / 55.0 };

template <auto Func_>
double stepRungeKutta45( double& inX, double& inY, double& inDif,
                         double inRelTol );
template <auto Func_>
std::vector<double> stepSIMLRungeKutta45( double& inX, std::vector<double>& inY,
                                          double& inDif, double inRelTol );

template <auto Func_>
double stepRungeKutta45( double& inX, double& inY, double& inDif,
                         double inRelTol )
{
    double lK[6];
    for ( int i = 0; i < 6; ++i )
    {
        double lTmpX = inX;
        lTmpX += inDif * gCoeffKX[i];
        double lTmpY = inY;
        for ( int j = 0; j < i; ++j ) { lTmpY += lK[j] * gCoeffKY[i][j]; }
        lK[i] = inDif * Func_( lTmpX, lTmpY );
    }

    double lY4 = inY;
    for ( int i = 0; i < 5; ++i ) { lY4 += lK[i] * gCoeffY4[i]; }
    double lY5 = inY;
    for ( int i = 0; i < 6; ++i ) { lY5 += lK[i] * gCoeffY5[i]; }

    double lError = std::abs( ( lY5 - lY4 ) / inDif );
    double lDelta = pow( 0.5 * inRelTol / lError, 0.25 );
    if ( lDelta < 0.1 ) { lDelta = 0.1; }
    if ( lDelta > 4.0 ) { lDelta = 4.0; }

    if ( lDelta > 1.0 )
    {
        inX += inDif;
        inDif *= lDelta;
        inY = lY5;
        return inY;
    }
    else
    {
        inDif *= lDelta;
        return stepRungeKutta45<Func_>( inX, inY, inDif, inRelTol );
    }
}

template <auto Func_>
std::pair<std::vector<double>, std::vector<double>> solveRungeKutta45(
    double inInitX, double inInitY, double inEndX, double inRelTol,
    double inMaxDif )
{
    std::vector<double> lXs, lYs;
    double lX   = inInitX;
    double lY   = inInitY;
    double lDif = ( inEndX - inInitX ) * 0.01;

    lXs.push_back( lX );
    lYs.push_back( lY );

    while ( ( lX - inEndX ) * ( inInitX - inEndX ) > 0 )
    {
        lY   = stepRungeKutta45<Func_>( lX, lY, lDif, inRelTol );
        lDif = std::min( { std::abs( lDif ), std::abs( inMaxDif ),
                           std::abs( lX - inEndX ) } );
        if ( inEndX < inInitX ) { lDif = -lDif; }
        lXs.push_back( lX );
        lYs.push_back( lY );
    }

    return std::make_pair( lXs, lYs );
}

template <auto Func_>
std::vector<double> stepSIMLRungeKutta45( double& inX, std::vector<double>& inY,
                                          double& inDif, double inRelTol )
{
    std::vector<std::vector<double>> lK( 6 );
    for ( int i = 0; i < 6; ++i )
    {
        double lTmpX = inX;
        lTmpX += inDif * gCoeffKX[i];
        std::vector<double> lTmpY = inY;
        for ( int j = 0; j < i; ++j )
        {
            for ( std::size_t iDim = 0; iDim < inY.size(); ++iDim )
            {
                lTmpY.at( iDim ) += lK.at( j ).at( iDim ) * gCoeffKY[i][j];
            }
        }
        lK.at( i ) = Func_( lTmpX, lTmpY );
        for ( std::size_t iDim = 0; iDim < inY.size(); ++iDim )
        {
            lK.at( i ).at( iDim ) *= inDif;
        }
    }

    std::vector<double> lY4 = inY;
    for ( int i = 0; i < 5; ++i )
    {
        for ( std::size_t iDim = 0; iDim < inY.size(); ++iDim )
        {
            lY4.at( iDim ) += lK.at( i ).at( iDim ) * gCoeffY4[i];
        }
    }
    std::vector<double> lY5 = inY;
    for ( int i = 0; i < 6; ++i )
    {
        for ( std::size_t iDim = 0; iDim < inY.size(); ++iDim )
        {
            lY5.at( iDim ) += lK.at( i ).at( iDim ) * gCoeffY5[i];
        }
    }

    double lError = 0.0;

    for ( std::size_t iDim = 0; iDim < inY.size(); ++iDim )
    {
        lError = std::max(
            lError, std::abs( ( lY5.at( iDim ) - lY4.at( iDim ) ) / inDif ) );
    }
    double lDelta = pow( 0.5 * inRelTol / lError, 0.25 );
    if ( lDelta < 0.1 ) { lDelta = 0.1; }
    if ( lDelta > 4.0 ) { lDelta = 4.0; }

    if ( lDelta > 1.0 )
    {
        inX += inDif;
        inDif *= lDelta;
        inY = lY5;
        return inY;
    }
    else
    {
        inDif *= lDelta;
        return stepSIMLRungeKutta45<Func_>( inX, inY, inDif, inRelTol );
    }
}

template <auto Func_>
std::vector<std::vector<double>> solveSIMLRungeKutta45(
    double inInitX, std::vector<double> inInitY, double inEndX, double inRelTol,
    double inMaxDif )
{
    std::vector<std::vector<double>> lResults( inInitY.size() + 1 );
    double lX              = inInitX;
    std::vector<double> lY = inInitY;
    double lDif            = ( inEndX - inInitX ) * 0.01;

    lResults.at( 0 ).push_back( lX );
    for ( std::size_t iDim = 0; iDim < inInitY.size(); ++iDim )
    {
        lResults.at( iDim + 1 ).push_back( lY.at( iDim ) );
    }

    while ( ( lX - inEndX ) * ( inInitX - inEndX ) > 0 )
    {
        lY   = stepSIMLRungeKutta45<Func_>( lX, lY, lDif, inRelTol );
        lDif = std::min( { std::abs( lDif ), std::abs( inMaxDif ),
                           std::abs( lX - inEndX ) } );
        if ( inEndX < inInitX ) { lDif = -lDif; }
        lResults.at( 0 ).push_back( lX );
        for ( std::size_t iDim = 0; iDim < inInitY.size(); ++iDim )
        {
            lResults.at( iDim + 1 ).push_back( lY.at( iDim ) );
        }
    }

    return lResults;
}

template <auto Func_>
std::vector<double> tmpFunctionForSecondOrderRungeKutta45(
    double inX, std::vector<double> inY )
{
    std::vector<double> lResults( 2 );
    lResults.at( 0 ) = inY.at( 1 );
    lResults.at( 1 ) = Func_( inX, inY.at( 0 ), inY.at( 1 ) );
    return lResults;
}

template <auto Func_>
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
solveSecondOrderRungeKutta45( double inInitX, double inInitY, double inInitDY,
                              double inEndX, double inRelTol, double inMaxDif )
{
    auto lResults =
        solveSIMLRungeKutta45<tmpFunctionForSecondOrderRungeKutta45<Func_>>(
            inInitX, { inInitY, inInitDY }, inEndX, inRelTol, inMaxDif );
    return std::make_tuple( lResults.at( 0 ), lResults.at( 1 ),
                            lResults.at( 2 ) );
}

template <auto Func_>
std::vector<double> tmpFunctionForSIMLSecondOrderRungeKutta45(
    double inX, std::vector<double> inY )
{
    std::vector<double> lResults( inY.size() );
    std::size_t lHalfSize = inY.size() / 2;
    std::vector<double> lY( inY.begin(), inY.begin() + lHalfSize );
    std::vector<double> lDY( inY.begin() + lHalfSize, inY.end() );
    std::vector<double> lTmpRes = Func_( inX, lY, lDY );
    for ( std::size_t i = 0; i < lHalfSize; ++i )
    {
        lResults.at( i )             = inY.at( lHalfSize + i );
        lResults.at( lHalfSize + i ) = lTmpRes.at( i );
    }
    return lResults;
}

template <auto Func_>
std::vector<std::vector<double>> solveSIMLSecondOrderRungeKutta45(
    double inInitX, std::vector<double> inInitY, std::vector<double> inInitDY,
    double inEndX, double inRelTol, double inMaxDif )
{
    std::copy( inInitDY.begin(), inInitDY.end(),
               std::back_inserter( inInitY ) );
    return solveSIMLRungeKutta45<
        tmpFunctionForSIMLSecondOrderRungeKutta45<Func_>>(
        inInitX, inInitY, inEndX, inRelTol, inMaxDif );
}

}  // namespace ODE
}  // namespace Math

#endif