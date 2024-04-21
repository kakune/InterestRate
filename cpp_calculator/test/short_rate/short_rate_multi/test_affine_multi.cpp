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

ShortRate::MultiFactor::CIR2ppWithMarket G2ppBuild(
    std::size_t inNTerms, std::size_t inNPath, double inMaturity,
    Math::Vec inConvSH, Math::Vec inMean, Math::Vec inVol,
    const Process::MarketData::ZCB& inMarketZCB )
{
    ShortRate::MultiFactor::CIR2ppWithMarketBuilder lBuilder;
    auto lTerms = makeTerms( inNTerms, inMaturity );
    lBuilder.setNPath( inNPath );
    lBuilder.setDrift( inConvSH, inMean );
    lBuilder.setVol( inVol );
    lBuilder.setTerms( lTerms );
    lBuilder.setMarketZCB( inMarketZCB );
    lBuilder.setRandom(
        std::make_unique<Process::RandomVec::StdBrownAntithetic>( 2 ) );
    return lBuilder.build();
}

double testCIR2ppConsistencyZCB( std::size_t inNTerms, std::size_t inNPath,
                                 double inMaturity, Math::Vec inConvSH,
                                 Math::Vec inMean, Math::Vec inVol,
                                 double inMeanDriftZCB, double inVolDriftZCB,
                                 std::size_t inSeed = 0 )
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
    auto lObj = G2ppBuild( inNTerms, inNPath, inMaturity, inConvSH, inMean,
                           inVol, lMarketZCB );
    Process::MarketData::ZCB lG2pp = lObj.createSpotRates().getZCB();

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

TEST( ShortRateMultiConstantGaussTest, CIR2ppZCBConsistency )
{
    EXPECT_NEAR(
        0.0,
        testCIR2ppConsistencyZCB( 100, 10000, 1.0, { 0.1, 0.5 }, { 0.03, 0.1 },
                                  { 0.02, 0.04 }, 0.1, 0.02 ),
        0.001 );
    EXPECT_NEAR(
        0.0,
        testCIR2ppConsistencyZCB( 100, 10000, 1.0, { 0.1, 0.5 }, { 0.03, 0.1 },
                                  { 0.02, 0.04 }, 0.1, 0.03 ),
        0.001 );
    EXPECT_NEAR(
        0.0,
        testCIR2ppConsistencyZCB( 200, 20000, 10.0, { 0.1, 0.15 },
                                  { 0.2, 0.06 }, { 0.01, 0.01 }, 0.2, 0.03 ),
        0.005 );
}
