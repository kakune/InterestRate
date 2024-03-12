#include <gtest/gtest.h>

#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "process/market.hpp"
#include "process/random.hpp"
#include "process/short_rate.hpp"

double testDifAnalytical( std::size_t inNTerms, std::size_t inNPath,
                          double inMaturity, double inInitialRate,
                          double inVol )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 1, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 1; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms  = std::make_shared<std::vector<double> >( lTerms );
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        inNPath, lsTerms );

    Process::ShortRate::HoLeeBuilder lBuilder;
    lBuilder.setTerms( lsTerms );
    lBuilder.setNPath( inNPath );
    lBuilder.setVol( inVol );
    lBuilder.setInitSpotRate( inInitialRate );
    lBuilder.setRandom( std::move( luRandom ) );

    Process::ShortRate::HoLee lHoLee = lBuilder.build();
    lHoLee.build();

    double lResult = 0.0;
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        double lNumerical  = lHoLee.priceZCB( 0, iTerm );
        double lAnalytical = lHoLee.analyticalPriceZCB( 0, iTerm );
        // std::cout << lNumerical << "," << lAnalytical << std::endl;
        lResult = std::max(
            lResult, std::abs( lNumerical - lAnalytical ) / lAnalytical );
    }
    return lResult;
}

double testConsistencyZCB( std::size_t inNTerms, std::size_t inNPath,
                           double inMaturity, double inVol,
                           double inMeanDriftZCB, double inVolDriftZCB,
                           std::size_t inSeed = 0 )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 2, 0 );
    std::vector<double> lZCB( inNTerms + 2, 1.0 );

    std::mt19937_64 lEngine( inSeed );
    std::uniform_real_distribution<double> lRandomGen(
        ( inMeanDriftZCB - inVolDriftZCB ) / double( inNTerms + 2 ),
        ( inMeanDriftZCB + inVolDriftZCB ) / double( inNTerms + 2 ) );

    for ( std::size_t iTerm = 1; iTerm < inNTerms + 2; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
        lZCB[iTerm]   = lZCB[iTerm - 1] - lRandomGen( lEngine );
    }
    auto lsTerms  = std::make_shared<std::vector<double> >( lTerms );
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        inNPath, lsTerms );

    Process::Market::Data lMarketData( lsTerms );
    lMarketData.setZCB( lZCB );

    Process::ShortRate::HoLeeBuilder lBuilder;
    lBuilder.setTerms( lsTerms );
    lBuilder.setNPath( inNPath );
    lBuilder.setVol( inVol );
    lBuilder.setMarketData(
        std::make_shared<Process::Market::Data>( lMarketData ) );
    lBuilder.setRandom( std::move( luRandom ) );

    Process::ShortRate::HoLee lHoLee = lBuilder.build();
    lHoLee.build();

    double lResult = 0.0;

    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        // std::cout << lZCB[iTerm] << "," << lHoLee.priceZCB( 0, iTerm )
        //           << std::endl;
        lResult = std::max(
            lResult, std::abs( ( lZCB[iTerm] - lHoLee.priceZCB( 0, iTerm ) ) /
                               lZCB[iTerm] ) );
    }
    return lResult;
}

TEST( ShortRateHoLeeTest, DifAnalytical )
{
    EXPECT_NEAR( 0.0, testDifAnalytical( 100, 10000, 1.0, 0.1, 0.05 ), 0.001 );
}

TEST( ShortRateHoLeeTest, ConsistencyZCB )
{
    EXPECT_NEAR(
        0.0, testConsistencyZCB( 100, 100000, 1.0, 0.05, 0.1, 0.05, 0 ), 0.01 );
    EXPECT_NEAR(
        0.0, testConsistencyZCB( 100, 100000, 1.0, 0.15, 0.1, 0.05, 0 ), 0.01 );
}