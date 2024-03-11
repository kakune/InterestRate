#include <gtest/gtest.h>

#include <iostream>
#include <random>
#include <vector>

#include "process/market.hpp"

double testConsistencyZCB( std::size_t inNTerms, double inMaturity,
                           double inMean, double inVol, std::size_t inSeed = 0 )
{
    double lDt = inMaturity / double( inNTerms );
    std::vector<double> lTerms( inNTerms, 0 );
    std::vector<double> lZCB( inNTerms, 1.0 );

    std::mt19937_64 lEngine( inSeed );
    std::uniform_real_distribution<double> lRandomGen(
        ( inMean - inVol ) / double( inNTerms ),
        ( inMean + inVol ) / double( inNTerms ) );

    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
        lZCB[iTerm]   = lZCB[iTerm - 1] - lRandomGen( lEngine );
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::Market::Data lObj( lsTerms );
    Process::Market::Data lObj2( lsTerms );

    lObj.setZCB( lZCB );
    std::vector<double> lFR( inNTerms );
    for ( std::size_t iTerm = 0; iTerm < inNTerms; ++iTerm )
    {
        lFR.at( iTerm ) = lObj.mInterpInstantaneousForwardRate( lTerms[iTerm] );
    }
    lObj2.setInstantaneousForwardRate( lFR );

    double lResult = 0.0;
    for ( std::size_t iTerm = 0; iTerm < inNTerms; ++iTerm )
    {
        // std::cout << lZCB[iTerm] << "," << lObj2.mInterpZCB( lTerms[iTerm] )
        //           << std::endl;
        lResult = std::max(
            0.0, std::abs( ( lZCB[iTerm] - lObj2.mInterpZCB( lTerms[iTerm] ) ) /
                           lZCB[iTerm] ) );
    }
    return lResult;
}

TEST( MarketTest, ConsistencyZCB )
{
    EXPECT_NEAR( 0.0, testConsistencyZCB( 100, 1.0, 0.1, 0.05, 0 ), 0.001 );
    EXPECT_NEAR( 0.0, testConsistencyZCB( 100, 1.0, 0.1, 0.02, 0 ), 0.001 );
}