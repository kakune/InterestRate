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

double testDifVasicek( std::size_t inNTerms, std::size_t inNPath,
                       double inMaturity, double inInitRate, double inVol,
                       double inKappa, double inMean )
{
    auto lTerms = makeTerms( inNTerms, inMaturity );
    auto luRandomVasicek =
        std::make_unique<Process::Random::StdBrownAntithetic>();
    auto luRandomGSR = std::make_unique<Process::Random::StdBrownAntithetic>();
    std::vector<double> lVols( inNTerms + 2, inVol );
    std::vector<double> lKappas( inNTerms + 2, inKappa );
    std::vector<double> lMeans( inNTerms + 2, inMean );

    ShortRate::OneFactor::VasicekBuilder lBuilderVasicek;
    lBuilderVasicek.setTerms( lTerms );
    lBuilderVasicek.setNPath( inNPath );
    lBuilderVasicek.setRandom( std::move( luRandomVasicek ) );
    lBuilderVasicek.setInitSpotRate( inInitRate );
    lBuilderVasicek.setVol( inVol );
    lBuilderVasicek.setKappa( inKappa );
    lBuilderVasicek.setMean( inMean );

    ShortRate::OneFactor::GSRBuilder lBuilderGSR;
    lBuilderGSR.setTerms( lTerms );
    lBuilderGSR.setNPath( inNPath );
    lBuilderGSR.setRandom( std::move( luRandomGSR ) );
    lBuilderGSR.setInitSpotRate( inInitRate );
    lBuilderGSR.setInterpVol( lVols );
    lBuilderGSR.setInterpKappa( lKappas );
    lBuilderGSR.setInterpMean( lMeans );

    ShortRate::OneFactor::Vasicek lVasicek = lBuilderVasicek.build();
    ShortRate::OneFactor::GSR lGSR         = lBuilderGSR.build();

    Process::MarketData::ZCB lVasicekZCB = lVasicek.createSpotRates().getZCB();
    Process::MarketData::ZCB lGSRZCB     = lGSR.createSpotRates().getZCB();

    double lResult = 0.0;
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lResult = std::max(
            { lResult, std::abs( lVasicekZCB[iTerm] - lGSRZCB[iTerm] ) /
                           lGSRZCB[iTerm] } );
    }
    return lResult;
}

double testDifVasicekWithMarket( std::size_t inNTerms, std::size_t inNPath,
                                 double inMaturity, double inVol,
                                 double inKappa, double inMeanDriftZCB,
                                 double inVolDriftZCB, std::size_t inSeed = 0 )
{
    auto lTerms = makeTerms( inNTerms, inMaturity );
    auto luRandomVasicek =
        std::make_unique<Process::Random::StdBrownAntithetic>();
    auto luRandomGSR = std::make_unique<Process::Random::StdBrownAntithetic>();
    std::vector<double> lVols( inNTerms + 2, inVol );
    std::vector<double> lKappas( inNTerms + 2, inKappa );
    std::vector<double> lZCB( inNTerms + 2, 1.0 );
    std::mt19937_64 lEngine( inSeed );
    std::uniform_real_distribution<double> lRandomGen(
        ( inMeanDriftZCB - inVolDriftZCB ) / double( inNTerms ),
        ( inMeanDriftZCB + inVolDriftZCB ) / double( inNTerms ) );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 2; ++iTerm )
    {
        lZCB[iTerm] = lZCB[iTerm - 1] - lRandomGen( lEngine );
    }
    Process::MarketData::ZCB lMarketZCB( lTerms, lZCB );

    ShortRate::OneFactor::VasicekWithMarketBuilder lBuilderVasicek;
    lBuilderVasicek.setTerms( lTerms );
    lBuilderVasicek.setNPath( inNPath );
    lBuilderVasicek.setRandom( std::move( luRandomVasicek ) );
    lBuilderVasicek.setVol( inVol );
    lBuilderVasicek.setKappa( inKappa );
    lBuilderVasicek.setMarketZCB( lMarketZCB );

    ShortRate::OneFactor::GSRWithMarketBuilder lBuilderGSR;
    lBuilderGSR.setTerms( lTerms );
    lBuilderGSR.setNPath( inNPath );
    lBuilderGSR.setRandom( std::move( luRandomGSR ) );
    lBuilderGSR.setInterpVol( lVols );
    lBuilderGSR.setInterpKappa( lKappas );
    lBuilderGSR.setMarketZCB( lMarketZCB );

    ShortRate::OneFactor::VasicekWithMarket lVasicek = lBuilderVasicek.build();
    ShortRate::OneFactor::GSRWithMarket lGSR         = lBuilderGSR.build();

    Process::MarketData::ZCB lVasicekZCB = lVasicek.createSpotRates().getZCB();
    Process::MarketData::ZCB lGSRZCB     = lGSR.createSpotRates().getZCB();

    double lResult = 0.0;
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lResult = std::max(
            { lResult, std::abs( lVasicekZCB[iTerm] - lGSRZCB[iTerm] ) /
                           lGSRZCB[iTerm] } );
    }
    return lResult;
}

TEST( ShortRateGSRTest, DifVasicek )
{
    EXPECT_NEAR( 0.0, testDifVasicek( 40, 10000, 1.0, 0.1, 0.1, 0.05, 0.05 ),
                 0.001 );
    EXPECT_NEAR( 0.0, testDifVasicek( 40, 10000, 1.0, 0.1, 0.2, 0.5, 0.2 ),
                 0.001 );
    EXPECT_NEAR( 0.0, testDifVasicek( 40, 10000, 1.0, 0.1, 0.3, 0.1, 0.4 ),
                 0.001 );
}

TEST( ShortRateGSRTest, DifVasicekWithMarket )
{
    EXPECT_NEAR(
        0.0,
        testDifVasicekWithMarket( 40, 100000, 1.0, 0.05, 0.1, 0.1, 0.005, 0 ),
        0.001 );
    EXPECT_NEAR(
        0.0,
        testDifVasicekWithMarket( 40, 100000, 1.0, 0.15, 0.1, 0.1, 0.005, 10 ),
        0.001 );
}