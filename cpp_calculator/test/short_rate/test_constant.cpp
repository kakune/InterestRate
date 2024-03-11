#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "process/short_rate.hpp"

double testPriceZCB( std::size_t inNTerms, std::size_t inNPath,
                     double inMaturity, double inRate )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 1, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 1; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::ShortRate::ConstantRate lObj( lsTerms, nullptr, inRate );
    lObj.build();
    return lObj.priceZCB( 0.0, inMaturity );
}

double testForwardRate( std::size_t inNTerms, std::size_t inNPath,
                        double inMaturity, double inRate, double inStartTime,
                        double inTerminalTime )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 1, 0.0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 1; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::ShortRate::ConstantRate lObj( lsTerms, nullptr, inRate );
    lObj.build();
    return lObj.forwardRate( inStartTime, inTerminalTime );
}

double testInstantaneousForwardRate( std::size_t inNTerms, std::size_t inNPath,
                                     double inMaturity, double inRate,
                                     double inFRTime )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 1, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 1; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::ShortRate::ConstantRate lObj( lsTerms, nullptr, inRate );
    lObj.build();
    return lObj.instantaneousForwardRate( inFRTime );
}

TEST( ShortRateConstantTest, PriceZCB )
{
    EXPECT_NEAR( std::exp( -0.1 ), testPriceZCB( 10, 10, 1.0, 0.1 ), 0.001 );
    EXPECT_NEAR( std::exp( -0.2 ), testPriceZCB( 10, 10, 1.0, 0.2 ), 0.001 );
    EXPECT_NEAR( std::exp( -6.0 ), testPriceZCB( 10, 10, 20.0, 0.3 ), 0.001 );
}

TEST( ShortRateConstantTest, ForwardRate )
{
    EXPECT_NEAR( 0.1, testForwardRate( 10, 10, 1.0, 0.1, 0.2, 0.5 ), 0.001 );
    EXPECT_NEAR( 0.2, testForwardRate( 10, 10, 1.0, 0.2, 0.3, 0.7 ), 0.001 );
    EXPECT_NEAR( 0.3, testForwardRate( 10, 10, 10.0, 0.3, 0.6, 10.0 ), 0.001 );
}

TEST( ShortRateConstantTest, InstantaneousForwardRate )
{
    EXPECT_NEAR( 0.1, testInstantaneousForwardRate( 100, 10, 1.0, 0.1, 0.5 ),
                 0.001 );
    EXPECT_NEAR( 0.2, testInstantaneousForwardRate( 100, 10, 1.0, 0.2, 0.3 ),
                 0.001 );
    EXPECT_NEAR( 0.3, testInstantaneousForwardRate( 100, 10, 10.0, 0.3, 6.7 ),
                 0.001 );
}