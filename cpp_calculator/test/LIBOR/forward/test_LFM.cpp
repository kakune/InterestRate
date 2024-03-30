#include <gtest/gtest.h>

#include "LIBOR/forward.hpp"

Process::MarketData::Terms makeTerms( std::size_t inNTerms, double inMaturity )
{
    double lDt = inMaturity * 1.0 / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    return Process::MarketData::Terms( lTerms );
}

LIBOR::Forward::ConstantLFM LFMBuild( std::size_t inNTerms, std::size_t inNPath,
                                      double inMaturity,
                                      std::vector<std::size_t> inIndTenor,
                                      Math::Vec inInitFRs, Math::Vec inVol,
                                      Math::Mat inCorr )
{
    LIBOR::Forward::ConstantLFMBuilder lBuilder;
    auto lTerms = makeTerms( inNTerms, inMaturity );
    lBuilder.setNPath( inNPath );
    lBuilder.setIndTenor( inIndTenor );
    lBuilder.setTerms( lTerms );
    lBuilder.setVols( inVol );
    lBuilder.setInitFRs( inInitFRs );
    lBuilder.setRandom(
        std::make_unique<Process::RandomVec::PathBrownAntithetic>( inNPath,
                                                                   lTerms, 2 ),
        inCorr );
    return lBuilder.build();
}

double testConstant( std::size_t inNTerms, std::size_t inNPath,
                     double inMaturity, std::vector<std::size_t> inIndTenor,
                     Math::Vec inInitFRs, Math::Vec inVol, Math::Mat inCorr )
{
    auto lObj = LFMBuild( inNTerms, inNPath, inMaturity, inIndTenor, inInitFRs,
                          inVol, inCorr );
    auto lFR  = lObj.createForwardRates();
    return 0.0;
}

TEST( ShortRateConstantTest, PriceZCB )
{
    std::vector<std::size_t> lIndTenor{ 100, 200, 300 };
    Math::Vec lInitFRs{ 0.1, 0.2 };
    Math::Vec lVol{ 0.9, 0.5 };
    Math::Mat lCorr{ { 1.0, 0.0 }, { 0.0, 1.0 } };
    EXPECT_NEAR(
        0.0, testConstant( 300, 10000, 1.0, lIndTenor, lInitFRs, lVol, lCorr ),
        0.0001 );
}
