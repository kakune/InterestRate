#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "process/short_rate_MC.hpp"

double testConstantPriceZCB( std::size_t inNTerms, std::size_t inNPath,
                             double inMaturity, double inRate )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTermsVec( inNTerms + 1, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 1; ++iTerm )
    {
        lTermsVec[iTerm] = lTermsVec[iTerm - 1] + lDt;
    }
    auto lTerms          = Process::MarketData::Terms( lTermsVec );
    Math::Vec lInitState = { inRate, 1.0 };
    Process::ShortRateMCMulti::ConstantRate lObj( lTerms, lInitState );
    return Process::MarketData::ZCB( lObj.calcSpotRates() )( 0.0, inMaturity );
}

double testConstantForwardRate( std::size_t inNTerms, std::size_t inNPath,
                                double inMaturity, double inRate,
                                double inStartTime, double inTerminalTime )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTermsVec( inNTerms + 1, 0.0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 1; ++iTerm )
    {
        lTermsVec[iTerm] = lTermsVec[iTerm - 1] + lDt;
    }
    auto lTerms          = Process::MarketData::Terms( lTermsVec );
    Math::Vec lInitState = { inRate, 1.0 };
    Process::ShortRateMCMulti::ConstantRate lObj( lTerms, lInitState );
    return Process::MarketData::ZCB( lObj.calcSpotRates() )
        .forwardRate( inStartTime, inTerminalTime );
}

double testConstantInstantaneousForwardRate( std::size_t inNTerms,
                                             std::size_t inNPath,
                                             double inMaturity, double inRate,
                                             double inFRTime )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTermsVec( inNTerms + 1, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 1; ++iTerm )
    {
        lTermsVec[iTerm] = lTermsVec[iTerm - 1] + lDt;
    }
    auto lTerms          = Process::MarketData::Terms( lTermsVec );
    Math::Vec lInitState = { inRate, 1.0 };
    Process::ShortRateMCMulti::ConstantRate lObj( lTerms, lInitState );
    return Process::MarketData::ZCB( lObj.calcSpotRates() )
        .instantaneousForwardRate( inFRTime );
}

TEST( ShortRateConstantTest, PriceZCB )
{
    EXPECT_NEAR( std::exp( -0.1 ), testConstantPriceZCB( 10, 10, 1.0, 0.1 ),
                 0.001 );
    EXPECT_NEAR( std::exp( -0.2 ), testConstantPriceZCB( 10, 10, 1.0, 0.2 ),
                 0.001 );
    EXPECT_NEAR( std::exp( -6.0 ), testConstantPriceZCB( 10, 10, 20.0, 0.3 ),
                 0.001 );
}

TEST( ShortRateConstantTest, ForwardRate )
{
    EXPECT_NEAR( 0.1, testConstantForwardRate( 10, 10, 1.0, 0.1, 0.2, 0.5 ),
                 0.001 );
    EXPECT_NEAR( 0.2, testConstantForwardRate( 10, 10, 1.0, 0.2, 0.3, 0.7 ),
                 0.001 );
    EXPECT_NEAR( 0.3, testConstantForwardRate( 10, 10, 10.0, 0.3, 0.6, 10.0 ),
                 0.001 );
}

TEST( ShortRateConstantTest, InstantaneousForwardRate )
{
    EXPECT_NEAR( 0.1,
                 testConstantInstantaneousForwardRate( 100, 10, 1.0, 0.1, 0.5 ),
                 0.001 );
    EXPECT_NEAR( 0.2,
                 testConstantInstantaneousForwardRate( 100, 10, 1.0, 0.2, 0.3 ),
                 0.001 );
    EXPECT_NEAR(
        0.3, testConstantInstantaneousForwardRate( 100, 10, 10.0, 0.3, 6.7 ),
        0.001 );
}