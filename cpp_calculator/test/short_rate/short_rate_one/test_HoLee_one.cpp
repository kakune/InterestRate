#include <gtest/gtest.h>

#include "short_rate/one-factor.hpp"

Process::MarketData::Terms makeTerms( std::size_t inNTerms, double inMaturity )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 2, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 2; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    return Process::MarketData::Terms( lTerms );
}

double testDifAnalytical( std::size_t inNTerms, std::size_t inNPath,
                          double inMaturity, double inInitialRate,
                          double inVol )
{
    auto lTerms   = makeTerms( inNTerms, inMaturity );
    auto luRandom = std::make_unique<Process::Random::StdBrownAntithetic>();

    ShortRate::OneFactor::HoLeeBuilder lBuilder;
    lBuilder.setTerms( lTerms );
    lBuilder.setNPath( inNPath );
    lBuilder.setVol( inVol );
    lBuilder.setInitSpotRate( inInitialRate );
    lBuilder.setRandom( std::move( luRandom ) );

    ShortRate::OneFactor::HoLee lHoLee = lBuilder.build();
    Process::MarketData::ZCB lHoLeeZCB = lHoLee.createSpotRates().getZCB();

    double lResult = 0.0;
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        double lNumerical  = lHoLeeZCB( 0.0, lTerms[iTerm] );
        double lAnalytical = lHoLee.analyticalPriceZCB( 0.0, lTerms[iTerm] );
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
    std::vector<double> lZCB( inNTerms + 2, 1.0 );
    std::mt19937_64 lEngine( inSeed );
    std::uniform_real_distribution<double> lRandomGen(
        ( inMeanDriftZCB - inVolDriftZCB ) / double( inNTerms ),
        ( inMeanDriftZCB + inVolDriftZCB ) / double( inNTerms ) );

    for ( std::size_t iTerm = 1; iTerm < inNTerms + 2; ++iTerm )
    {
        lZCB[iTerm] = lZCB[iTerm - 1] - lRandomGen( lEngine );
    }

    auto lTerms   = makeTerms( inNTerms, inMaturity );
    auto luRandom = std::make_unique<Process::Random::StdBrownAntithetic>();

    Process::MarketData::ZCB lMarketZCB( lTerms, lZCB );

    ShortRate::OneFactor::HoLeeWithMarketBuilder lBuilder;
    lBuilder.setTerms( lTerms );
    lBuilder.setNPath( inNPath );
    lBuilder.setVol( inVol );
    lBuilder.setMarketZCB( lMarketZCB );
    lBuilder.setRandom( std::move( luRandom ) );

    ShortRate::OneFactor::HoLeeWithMarket lHoLee = lBuilder.build();
    Process::MarketData::ZCB lHoLeeZCB = lHoLee.createSpotRates().getZCB();

    double lResult = 0.0;

    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        // std::cout << lZCB[iTerm] << "," << lHoLeeZCB[iTerm] << std::endl;
        lResult = std::max(
            lResult,
            std::abs( ( lZCB[iTerm] - lHoLeeZCB[iTerm] ) / lZCB[iTerm] ) );
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