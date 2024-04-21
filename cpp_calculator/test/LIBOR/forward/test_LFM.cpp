#include <gtest/gtest.h>

#include "LIBOR/forward.hpp"
#include "analytical/Black76.hpp"

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

// Following requirement cannot be interpreted by icpx. g++ can compile.
// template <template <LIBOR::Forward::C_VolGen> class StepCalculator_>
//     requires LIBOR::Forward::C_StepCalc<
//         StepCalculator_<LIBOR::Forward::VolGen::Constant>>
template <template <LIBOR::Forward::C_VolGen> class StepCalculator_>
double testImpVolConstantVol( std::size_t inNPath, double inMaturity,
                              std::vector<std::size_t> inIndTenor,
                              Math::Vec inInitFR, Math::Vec inVol,
                              Math::Mat inCorr, Math::Vec inCorrectImpVol )
{
    auto lTerms = makeTerms( inIndTenor.back() + 1, inMaturity );
    auto lTenor = Process::MarketData::Tenor( lTerms, inIndTenor );

    LIBOR::Forward::VolGen::Constant lVolGen( inVol );
    StepCalculator_ lStep( lTerms, lTenor, inCorr, lVolGen );
    auto lFR =
        LIBOR::Forward::Factory( inNPath, lTerms, lTenor, inInitFR, lStep )
            .template createForwardRates<
                Process::RandomVec::StdBrownAntithetic>();
    for ( auto& fr : lFR[0] ) { fr.print(); }
    Math::Vec lImpVolByCaplet( inCorrectImpVol ),
        lImpVolByFloorlet( inCorrectImpVol );
    for ( std::size_t i = 1; i < inIndTenor.size() - 1; ++i )
    {
        lImpVolByCaplet( i )   = lFR.calcBlackImpVol( inInitFR( i ), i, true );
        lImpVolByFloorlet( i ) = lFR.calcBlackImpVol( inInitFR( i ), i, false );
    }
    std::cout << "implied volatility by caplet   : ", lImpVolByCaplet.print();
    std::cout << "implied volatility by floorlet : ", lImpVolByFloorlet.print();
    return std::max( abs( lImpVolByCaplet - inCorrectImpVol ).max(),
                     abs( lImpVolByFloorlet - inCorrectImpVol ).max() );
}

TEST( ShortRateConstantTest, PriceZCB )
{
    std::vector<std::size_t> lIndTenor{ 0, 1, 3, 4 };
    Math::Vec lInitFR{ 0.2, 0.2, 0.2 };
    Math::Vec lVol{ 0.01, 0.02, 0.05 };
    Math::Mat lCorr    = Math::unitMat( 3, 1.0 );
    std::size_t lNPath = 400000;
    double lMaturity   = 1.0;

    EXPECT_NEAR(
        0.0,
        testImpVolConstantVol<LIBOR::Forward::StepCalc::LogNormalTerminalMeas>(
            lNPath, lMaturity, lIndTenor, lInitFR, lVol, lCorr, lVol ),
        0.001 );
    EXPECT_NEAR( 0.0,
                 testImpVolConstantVol<
                     LIBOR::Forward::StepCalc::LogNormalTerminalMeasWithLog>(
                     lNPath, lMaturity, lIndTenor, lInitFR, lVol, lCorr, lVol ),
                 0.001 );
    EXPECT_NEAR(
        0.0,
        testImpVolConstantVol<LIBOR::Forward::StepCalc::LogNormalSpotMeas>(
            lNPath, lMaturity, lIndTenor, lInitFR, lVol, lCorr, lVol ),
        0.001 );
}
