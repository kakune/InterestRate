#include <gtest/gtest.h>

#include "short_rate/multi-factor.hpp"

Process::MarketData::Terms makeTerms( std::size_t inNTerms, double inMaturity )
{
    double lDt = inMaturity * 1.05 / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    return Process::MarketData::Terms( lTerms );
}

ShortRate::MultiFactor::ConstantGauss rateBuild(
    std::size_t inNTerms, std::size_t inNPath, double inMaturity,
    Math::Vec inInitState, Math::Vec inDriftCoeff, Math::Mat inVolCoeff )
{
    ShortRate::MultiFactor::ConstantGaussBuilder lBuilder;
    auto lTerms = makeTerms( inNTerms, inMaturity );
    lBuilder.setInitState( inInitState );
    lBuilder.setNPath( inNPath );
    lBuilder.setDrift( inDriftCoeff );
    lBuilder.setVol( inVolCoeff );
    lBuilder.setTerms( lTerms );
    lBuilder.setRandom(
        std::make_unique<Process::RandomVec::PathBrownAntithetic>(
            inNPath, lTerms, inVolCoeff.sizeCol() ) );
    return lBuilder.build();
}

ShortRate::MultiFactor::G2ppWithMarket G2ppBuild(
    std::size_t inNTerms, std::size_t inNPath, double inMaturity,
    Math::Vec inDriftCoeff, Math::Mat inVolCoeff,
    const Process::MarketData::ZCB& inMarketZCB )
{
    ShortRate::MultiFactor::G2ppWithMarketBuilder lBuilder;
    auto lTerms = makeTerms( inNTerms, inMaturity );
    lBuilder.setNPath( inNPath );
    lBuilder.setDrift( inDriftCoeff );
    lBuilder.setVol( inVolCoeff );
    lBuilder.setTerms( lTerms );
    lBuilder.setMarketZCB( inMarketZCB );
    lBuilder.setRandom(
        std::make_unique<Process::RandomVec::PathBrownAntithetic>(
            inNPath, lTerms, 2 ) );
    return lBuilder.build();
}

double testConstantPriceZCB( std::size_t inNTerms, std::size_t inNPath,
                             double inMaturity, Math::Vec inInitState,
                             Math::Vec inDriftCoeff, Math::Mat inVolCoeff )
{
    auto lObj = rateBuild( inNTerms, inNPath, inMaturity, inInitState,
                           inDriftCoeff, inVolCoeff );
    return lObj.createSpotRates().createZCB()( 0.0, inMaturity );
}

double testConstantPriceZCBAnalytical( double inMaturity, Math::Vec inInitState,
                                       Math::Vec inDriftCoeff,
                                       Math::Mat inVolCoeff )
{
    auto lObj =
        rateBuild( 1, 1, inMaturity, inInitState, inDriftCoeff, inVolCoeff );
    return lObj.analyticalPriceZCB( inMaturity );
}

double testG2ppConsistencyZCB( std::size_t inNTerms, std::size_t inNPath,
                               double inMaturity, Math::Vec inDriftCoeff,
                               Math::Mat inVolCoeff, double inMeanDriftZCB,
                               double inVolDriftZCB, std::size_t inSeed = 0 )
{
    std::vector<double> lZCB( inNTerms, 1.0 );
    std::mt19937_64 lEngine( inSeed );
    std::uniform_real_distribution<double> lRandomGen(
        ( inMeanDriftZCB - inVolDriftZCB ) / double( inNTerms ),
        ( inMeanDriftZCB + inVolDriftZCB ) / double( inNTerms ) );

    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lZCB[iTerm] = lZCB[iTerm - 1] - lRandomGen( lEngine );
    }

    auto lTerms = makeTerms( inNTerms, inMaturity );
    Process::MarketData::ZCB lMarketZCB( lTerms, lZCB );
    auto lObj = G2ppBuild( inNTerms, inNPath, inMaturity, inDriftCoeff,
                           inVolCoeff, lMarketZCB );
    Process::MarketData::ZCB lG2pp = lObj.createSpotRates().createZCB();

    double lResult = 0.0;

    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        // std::cout << lZCB[iTerm] << ", " << lG2pp[iTerm] << std::endl;
        lResult = std::max( lResult,
                            std::abs( ( lZCB[iTerm] - lG2pp( lTerms[iTerm] ) ) /
                                      lZCB[iTerm] ) );
    }
    return lResult;
}

TEST( ShortRateMultiConstantGaussTest, ZCB )
{
    EXPECT_NEAR(
        testConstantPriceZCB( 100, 10000, 1.0, { 0.1, 0.05 }, { -0.01, -0.005 },
                              { { 0.03, 0.0 }, { 0.02, 0.04 } } ),
        testConstantPriceZCBAnalytical( 1.0, { 0.1, 0.05 }, { -0.01, -0.005 },
                                        { { 0.03, 0.0 }, { 0.02, 0.04 } } ),
        0.002 );
    EXPECT_NEAR(
        testConstantPriceZCB( 100, 10000, 1.0, { 0.1, 0.2 }, { -0.005, -0.01 },
                              { { 0.08, 0.0 }, { 0.09, 0.05 } } ),
        testConstantPriceZCBAnalytical( 1.0, { 0.1, 0.2 }, { -0.005, -0.01 },
                                        { { 0.08, 0.0 }, { 0.09, 0.05 } } ),
        0.002 );
}

TEST( ShortRateMultiConstantGaussTest, G2ppZCBConsistency )
{
    EXPECT_NEAR(
        0.0,
        testG2ppConsistencyZCB( 100, 10000, 1.0, { -0.01, -0.005 },
                                { { 0.03, 0.0 }, { 0.02, 0.04 } }, 0.1, 0.02 ),
        0.001 );
    EXPECT_NEAR(
        0.0,
        testG2ppConsistencyZCB( 100, 10000, 1.0, { -0.01, 0.005 },

                                { { 0.03, 0.0 }, { 0.02, -0.04 } }, 0.1, 0.03 ),
        0.001 );
    EXPECT_NEAR( 0.0,
                 testG2ppConsistencyZCB( 200, 20000, 10.0, { -0.01, -0.015 },
                                         { { -0.02, 0.0 }, { 0.01, -0.01 } },
                                         0.2, 0.03 ),
                 0.005 );
}
