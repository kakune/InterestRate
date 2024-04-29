#include <gtest/gtest.h>

#include <vector>

#include "analytical/Black76.hpp"
#include "process/market_data.hpp"

static Process::MarketData::Terms makeTerms( std::size_t inNTerms,
                                             double inMaturity )
{
    double lDt = inMaturity * 1.0 / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    return Process::MarketData::Terms( lTerms );
}

static Process::MarketData::ZCB makeBlackZCB(
    const Process::MarketData::Terms& inTerms, double inRate )
{
    std::vector<double> lZCBVec( inTerms.size(), 1.0 );
    for ( std::size_t iTerm = 1; iTerm < inTerms.size(); ++iTerm )
    {
        lZCBVec[iTerm] =
            lZCBVec[iTerm - 1] *
            std::exp( -inRate * ( inTerms[iTerm] - inTerms[iTerm - 1] ) );
    }
    return Process::MarketData::ZCB( inTerms, lZCBVec );
}

static Process::MarketData::Caplets makeBlackCaplets(
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::ZCB& inZCB, double inVol,
    const std::vector<double>& inStrikes )
{
    std::vector<double> lZCBVec( inTerms.size() - 1 );
    for ( std::size_t iTerm = 0; iTerm < lZCBVec.size(); ++iTerm )
    {
        lZCBVec[iTerm] = inZCB( inTerms[iTerm + 1] );
    }
    Analytical::Black76::Model lBlack( *( inTerms.ptr() ) );
    lBlack.setInitZCB( lZCBVec );
    lBlack.setVol( inVol );

    std::vector<std::vector<double>> lCapletVec(
        inStrikes.size(), std::vector<double>( inTerms.size() - 1 ) );
    for ( std::size_t iStrike = 0; iStrike < inStrikes.size(); ++iStrike )
    {
        for ( std::size_t iTerm = 0; iTerm < inTerms.size() - 1; ++iTerm )
        {
            lCapletVec[iStrike][iTerm] =
                lBlack.priceCaplet( inStrikes[iStrike], iTerm );
        }
    }
    return Process::MarketData::Caplets( inTerms, inStrikes, lCapletVec,
                                         inZCB );
}

double testInstantaneousForwardRateOfBlackZCB( std::size_t inNTerms,
                                               double inMaturity,
                                               double inRate )
{
    auto lTerms = makeTerms( inNTerms, inMaturity );
    auto lZCB   = makeBlackZCB( lTerms, inRate );

    double lResult = 0.0;
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lResult = std::max(
            lResult, std::abs( lZCB.instantaneousForwardRate( lTerms[iTerm] ) -
                               inRate ) );
    }
    return lResult;
}

double testImpVolOfBlackCaplet( std::size_t inNTerms, double inMaturity,
                                double inRate, double inVol,
                                std::vector<double> inStrikes )
{
    auto lTerms   = makeTerms( inNTerms, inMaturity );
    auto lZCB     = makeBlackZCB( lTerms, inRate );
    auto lCaplets = makeBlackCaplets( lTerms, lZCB, inVol, inStrikes );

    double lResult = 0.0;
    for ( std::size_t iStrike = 0; iStrike < inStrikes.size(); ++iStrike )
    {
        for ( std::size_t iTerm = 1; iTerm < lTerms.size() - 1; ++iTerm )
        {
            lResult = std::max( { lResult, std::abs( lCaplets.impliedBlackVol(
                                                         iStrike, iTerm ) -
                                                     inVol ) } );
        }
    }
    return lResult;
}

TEST( MarketTest, BlackZCB )
{
    EXPECT_NEAR( 0.0, testInstantaneousForwardRateOfBlackZCB( 100, 1.0, 0.1 ),
                 0.001 );
}

TEST( MarketTest, BlackCaplets )
{
    EXPECT_NEAR( 0.0, testImpVolOfBlackCaplet( 100, 1.0, 0.1, 0.25, { 0.1 } ),
                 0.001 );
}