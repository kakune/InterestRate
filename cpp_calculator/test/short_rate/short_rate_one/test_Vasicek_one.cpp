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
                          double inMaturity, double inInitRate, double inVol,
                          double inKappa, double inMean )
{
    auto lTerms   = makeTerms( inNTerms, inMaturity );
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        inNPath, lTerms );

    ShortRate::OneFactor::VasicekBuilder lBuilder;
    lBuilder.setTerms( lTerms );
    lBuilder.setRandom( std::move( luRandom ) );
    lBuilder.setNPath( inNPath );
    lBuilder.setInitSpotRate( inInitRate );
    lBuilder.setVol( inVol );
    lBuilder.setKappa( inKappa );
    lBuilder.setMean( inMean );

    ShortRate::OneFactor::Vasicek lVasicek = lBuilder.build();
    Process::MarketData::ZCB lVasicekZCB( lVasicek.calcSpotRates() );

    double lResult = 0.0;
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        double lNumerical  = lVasicekZCB( 0.0, lTerms[iTerm] );
        double lAnalytical = lVasicek.analyticalPriceZCB( 0, iTerm );
        // std::cout << lNumerical << "," << lAnalytical << std::endl;
        lResult = std::max(
            lResult, std::abs( lNumerical - lAnalytical ) / lAnalytical );
    }
    return lResult;
}

double testConsistencyZCB( std::size_t inNTerms, std::size_t inNPath,
                           double inMaturity, double inVol, double inKappa,
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
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        inNPath, lTerms );

    Process::MarketData::ZCB lMarketZCB( lTerms, lZCB );

    ShortRate::OneFactor::VasicekWithMarketBuilder lBuilder;
    lBuilder.setTerms( lTerms );
    lBuilder.setNPath( inNPath );
    lBuilder.setVol( inVol );
    lBuilder.setKappa( inKappa );
    lBuilder.setMarketZCB( lMarketZCB );
    lBuilder.setRandom( std::move( luRandom ) );

    ShortRate::OneFactor::VasicekWithMarket lVasicek = lBuilder.build();
    Process::MarketData::SpotRates lSpot             = lVasicek.calcSpotRates();
    Process::MarketData::ZCB lVasicekZCB( lSpot, 3 );

    double lResult = 0.0;

    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        std::cout << lTerms[iTerm] << "," << lZCB[iTerm] << ","
                  << lVasicekZCB[iTerm] << std::endl;
        lResult = std::max(
            lResult,
            std::abs( ( lZCB[iTerm] - lVasicekZCB[iTerm] ) / lZCB[iTerm] ) );
    }
    return lResult;
}

TEST( ShortRateVasicekTest, DifAnalytical )
{
    EXPECT_NEAR( 0, testDifAnalytical( 40, 50000, 1.0, 0.1, 0.05, 0.05, 0.1 ),
                 0.01 );
}

TEST( ShortRateVasicekTest, ConsistencyZCB )
{
    EXPECT_NEAR(
        0.0, testConsistencyZCB( 100, 100000, 1.0, 0.05, 0.1, 0.1, 0.005, 0 ),
        0.01 );
    EXPECT_NEAR(
        0.0, testConsistencyZCB( 100, 100000, 1.0, 0.15, 0.1, 0.1, 0.005, 10 ),
        0.01 );
}
