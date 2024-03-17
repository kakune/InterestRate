#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "process/short_rate.hpp"

std::shared_ptr<std::vector<double> > makeTerms( std::size_t inNTerms,
                                                 double inMaturity )
{
    double lDt = inMaturity / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms + 5, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms + 5; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    return std::make_shared<std::vector<double> >( lTerms );
}

Process::ShortRate::ConstantAffine rateBuild( std::size_t inNTerms,
                                              std::size_t inNPath,
                                              double inMaturity, double inRate,
                                              double inLambda, double inEta,
                                              double inGamma, double inDelta )
{
    Process::ShortRate::ConstantAffineBuilder lBuilder;
    auto lsTerms = makeTerms( inNTerms, inMaturity );
    lBuilder.setInitSpotRate( inRate );
    lBuilder.setNPath( 100 );
    lBuilder.setDrift( inLambda, inEta );
    lBuilder.setVol( inGamma, inDelta );
    lBuilder.setTerms( lsTerms );
    lBuilder.setRandom( std::make_unique<Process::Random::PathBrownAntithetic>(
        100, lsTerms ) );
    return lBuilder.build();
}
double testConstantPriceZCB( std::size_t inNTerms, std::size_t inNPath,
                             double inMaturity, double inRate )
{
    auto lObj =
        rateBuild( inNTerms, inNPath, inMaturity, inRate, 0.0, 0.0, 0.0, 0.0 );
    lObj.build();
    return lObj.priceZCB( 0.0, inMaturity );
}
double testConstantForwardRate( std::size_t inNTerms, std::size_t inNPath,
                                double inMaturity, double inRate,
                                double inStartTime, double inTerminalTime )
{
    auto lObj =
        rateBuild( inNTerms, inNPath, inMaturity, inRate, 0.0, 0.0, 0.0, 0.0 );
    lObj.build();
    return lObj.forwardRate( inStartTime, inTerminalTime );
}
double testConstantInstantaneousForwardRate( std::size_t inNTerms,
                                             std::size_t inNPath,
                                             double inMaturity, double inRate,
                                             double inFRTime )
{
    auto lObj =
        rateBuild( inNTerms, inNPath, inMaturity, inRate, 0.0, 0.0, 0.0, 0.0 );
    lObj.build();
    return lObj.instantaneousForwardRate( inFRTime );
}

double testCIRPriceZCB( std::size_t inNTerms, std::size_t inNPath,
                        double inMaturity, double inRate, double inK,
                        double inMean, double inVol )
{
    auto lObj = rateBuild( inNTerms, inNPath, inMaturity, inRate, -inK,
                           inK * inMean, inVol * inVol, 0.0 );
    lObj.build();
    return lObj.priceZCB( 0.0, inMaturity );
}
double testCIRInstantaneousForwardRate( std::size_t inNTerms,
                                        std::size_t inNPath, double inMaturity,
                                        double inRate, double inK,
                                        double inMean, double inVol )
{
    auto lObj = rateBuild( inNTerms, inNPath, inMaturity, inRate, -inK,
                           inK * inMean, inVol * inVol, 0.0 );
    lObj.build();
    return lObj.instantaneousForwardRate( inMaturity );
}

double testCIRPriceZCBByAB( std::size_t inNTerms, double inTimeMesh,
                            double inMaturity, double inRate, double inK,
                            double inMean, double inVol )
{
    auto lObj = rateBuild( inNTerms, 1, inMaturity, inRate, -inK, inK * inMean,
                           inVol * inVol, 0.0 );
    lObj.buildAB( inTimeMesh );
    return lObj.priceZCBByAB( inMaturity );
}
double testCIRInstantaneousForwardRateByAB( std::size_t inNTerms,
                                            double inTimeMesh,
                                            double inMaturity, double inRate,
                                            double inK, double inMean,
                                            double inVol )
{
    auto lObj = rateBuild( inNTerms, 1, inMaturity, inRate, -inK, inK * inMean,
                           inVol * inVol, 0.0 );
    lObj.buildAB( inTimeMesh );
    return lObj.instantaneousForwardRateByAB( inMaturity );
}

TEST( ShortRateConstantAffineTest, Constant )
{
    EXPECT_NEAR( std::exp( -0.1 ), testConstantPriceZCB( 10, 10, 1.0, 0.1 ),
                 0.001 );
    EXPECT_NEAR( std::exp( -0.2 ), testConstantPriceZCB( 10, 10, 1.0, 0.2 ),
                 0.001 );
    EXPECT_NEAR( std::exp( -6.0 ), testConstantPriceZCB( 10, 10, 20.0, 0.3 ),
                 0.001 );
    EXPECT_NEAR( 0.1, testConstantForwardRate( 10, 10, 1.0, 0.1, 0.2, 0.5 ),
                 0.001 );
    EXPECT_NEAR( 0.2, testConstantForwardRate( 10, 10, 1.0, 0.2, 0.3, 0.7 ),
                 0.001 );
    EXPECT_NEAR( 0.3, testConstantForwardRate( 10, 10, 10.0, 0.3, 0.6, 10.0 ),
                 0.001 );
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

TEST( ShortRateConstantAffineTest, CIR )
{
    EXPECT_NEAR( 0.987862017045984,
                 testCIRPriceZCB( 50, 1000, 1.0, 0.01, 0.05, 0.1, 0.02 ),
                 0.001 );
    EXPECT_NEAR( 0.3685782372351295,
                 testCIRPriceZCB( 50, 1000, 10.0, 0.1, 0.2, 0.1, 0.02 ),
                 0.001 );
    EXPECT_NEAR( 0.913980503378077,
                 testCIRPriceZCB( 50, 1000, 0.1, 0.9, 0.0, 0.7, 0.6 ), 0.001 );
    EXPECT_NEAR( 0.0143871638049957,
                 testCIRInstantaneousForwardRate( 100, 10000, 1.0, 0.01, 0.05,
                                                  0.1, 0.02 ),
                 0.001 );
    EXPECT_NEAR( 0.0996280717912332,
                 testCIRInstantaneousForwardRate( 100, 10000, 10.0, 0.1, 0.2,
                                                  0.1, 0.02 ),
                 0.001 );
    EXPECT_NEAR(
        0.898381942018978,
        testCIRInstantaneousForwardRate( 100, 10000, 0.1, 0.9, 0.0, 0.7, 0.6 ),
        0.005 );

    EXPECT_NEAR( 0.987862017045984,
                 testCIRPriceZCBByAB( 50.0, 50.0, 1.0, 0.01, 0.05, 0.1, 0.02 ),
                 0.001 );
    EXPECT_NEAR( 0.3685782372351295,
                 testCIRPriceZCBByAB( 50.0, 50.0, 10.0, 0.1, 0.2, 0.1, 0.02 ),
                 0.001 );
    EXPECT_NEAR( 0.913980503378077,
                 testCIRPriceZCBByAB( 50.0, 50.0, 0.1, 0.9, 0.0, 0.7, 0.6 ),
                 0.001 );
    EXPECT_NEAR( 0.0143871638049957,
                 testCIRInstantaneousForwardRateByAB( 50.0, 100.0, 1.0, 0.01,
                                                      0.05, 0.1, 0.02 ),
                 0.003 );
    EXPECT_NEAR( 0.0996280717912332,
                 testCIRInstantaneousForwardRateByAB( 50.0, 100.0, 10.0, 0.1,
                                                      0.2, 0.1, 0.02 ),
                 0.001 );
    EXPECT_NEAR( 0.898381942018978,
                 testCIRInstantaneousForwardRateByAB( 50.0, 100.0, 0.1, 0.9,
                                                      0.0, 0.7, 0.6 ),
                 0.005 );
}
