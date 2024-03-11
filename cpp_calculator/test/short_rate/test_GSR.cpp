#include <gtest/gtest.h>

#include <iostream>
#include <memory>
#include <vector>

#include "process/market.hpp"
#include "process/random.hpp"
#include "process/short_rate.hpp"

double testDifVasicek( std::size_t inNTerms, std::size_t inNPath,
                       double inMaturity, double inInitRate, double inVol,
                       double inKappa, double inMean )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 2, 0.0 );
    std::vector<double> lVols( inNTerms + 2, inVol );
    std::vector<double> lKappas( inNTerms + 2, inKappa );
    std::vector<double> lMeans( inNTerms + 2, inMean );

    for ( std::size_t iTerm = 1; iTerm < inNTerms + 2; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );
    auto luRandomVasicek =
        std::make_unique<Process::Random::PathBrownAntithetic>( inNPath,
                                                                lsTerms );
    auto luRandomGSR = std::make_unique<Process::Random::PathBrownAntithetic>(
        inNPath, lsTerms );

    Process::ShortRate::VasicekBuilder lBuilderVasicek;
    lBuilderVasicek.setTerms( lsTerms );
    lBuilderVasicek.setNPath( inNPath );
    lBuilderVasicek.setRandom( std::move( luRandomVasicek ) );
    lBuilderVasicek.setInitSpotRate( inInitRate );
    lBuilderVasicek.setVol( inVol );
    lBuilderVasicek.setKappa( inKappa );
    lBuilderVasicek.setMean( inMean );

    Process::ShortRate::GSRBuilder lBuilderGSR;
    lBuilderGSR.setTerms( lsTerms );
    lBuilderGSR.setNPath( inNPath );
    lBuilderGSR.setRandom( std::move( luRandomGSR ) );
    lBuilderGSR.setInitSpotRate( inInitRate );
    lBuilderGSR.setInterpVol( lVols );
    lBuilderGSR.setInterpKappa( lKappas );
    lBuilderGSR.setInterpMean( lMeans );

    Process::ShortRate::Vasicek lVasicek = lBuilderVasicek.build();
    Process::ShortRate::GSR lGSR         = lBuilderGSR.build();

    lVasicek.build();
    lGSR.build();

    double lResult = 0.0;
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        double lVasicekZCB = lVasicek.priceZCB( 0, iTerm );
        double lGSRZCB     = lGSR.priceZCB( 0, iTerm );
        double lErrZCB     = std::abs( lVasicekZCB - lGSRZCB ) / lGSRZCB;
        double lVasicekFR  = lVasicek.forwardRate( 0, iTerm );
        double lGSRFR      = lVasicek.forwardRate( 0, iTerm );
        double lErrFR      = std::abs( lVasicekFR - lGSRFR ) / lGSRFR;
        lResult            = std::max( { lResult, lErrZCB, lErrFR } );
    }

    return lResult;
}

TEST( ShortRateGSRTest, DifVasicek )
{
    EXPECT_NEAR( 0.0, testDifVasicek( 40, 50000, 1.0, 0.1, 0.1, 0.05, 0.05 ),
                 0.001 );
    EXPECT_NEAR( 0.0, testDifVasicek( 40, 50000, 1.0, 0.1, 0.2, 0.5, 0.2 ),
                 0.001 );
    EXPECT_NEAR( 0.0, testDifVasicek( 40, 50000, 1.0, 0.1, 0.3, 0.1, 0.4 ),
                 0.001 );
}