#include <gtest/gtest.h>

#include <iostream>
#include <random>
#include <vector>

#include "process/market.hpp"

bool testConsistencyZCB( double inMean, double inVol, std::size_t inSeed )
{
    std::size_t lNTerms = 100;
    double lMaturity    = 1.0;
    double lDt          = lMaturity / double( lNTerms );
    std::vector<double> lTerms( lNTerms, 0 );
    std::vector<double> lZCB( lNTerms, 1.0 );

    std::mt19937_64 lEngine( inSeed );
    std::uniform_real_distribution<double> lRandomGen(
        ( inMean - inVol ) / double( lNTerms ),
        ( inMean + inVol ) / double( lNTerms ) );

    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
        lZCB[iTerm]   = lZCB[iTerm - 1] - lRandomGen( lEngine );
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::Market::Data lObj( lsTerms );
    lObj.setZCB( lZCB );
    std::vector<double> lFR( lNTerms );
    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        lFR.at( iTerm ) = lObj.mInterpInstantaneousForwardRate( lTerms[iTerm] );
    }
    Process::Market::Data lObj2( lsTerms );
    lObj2.setForwardRate( lFR );
    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        if ( std::abs( ( lZCB[iTerm] - lObj2.mInterpZCB( lTerms[iTerm] ) ) /
                       lZCB[iTerm] ) > 0.01 )
        {
            return false;
        }
    }
    return true;
}

TEST( MarketTest, ConsistencyZCB )
{
    EXPECT_TRUE( testConsistencyZCB( 0.1, 0.05, 0 ) );
    EXPECT_TRUE( testConsistencyZCB( 0.1, 0.02, 1 ) );
    EXPECT_TRUE( testConsistencyZCB( 0.1, 0.08, 2 ) );
}