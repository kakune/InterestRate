#include <gtest/gtest.h>

#include "LIBOR/forward.hpp"
#include "process/random.hpp"

Process::MarketData::Terms makeTerms( std::size_t inNTerms, double inMaturity )
{
    double lDt = inMaturity * 1.0 / double( inNTerms - 1 );
    std::vector<double> lTerms( inNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < inNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
        // std::cout << iTerm << ":" << lTerms[iTerm] << std::endl;
    }
    return Process::MarketData::Terms( lTerms );
}

// Following requirement cannot be interpreted by icpx. g++ can compile.
// template <template <LIBOR::Forward::C_VolGen> class StepCalculator_>
//     requires LIBOR::Forward::C_StepCalc<
//         StepCalculator_<LIBOR::Forward::VolGen::Constant>>
template <template <LIBOR::Forward::C_VolGen> class StepCalculator_>
double testImpVolSABR( std::size_t inNPath, double inMaturity,
                       std::vector<std::size_t> inIndTenor, Math::Mat inCorr,
                       Math::Vec inInitFR, Math::Vec inInitVol,
                       double inExponent, double inVolVol, Math::Vec inCorrSV )
{
    auto lTerms = makeTerms( inIndTenor.back() + 1, inMaturity );
    auto lTenor = Process::MarketData::Tenor( lTerms, inIndTenor );

    LIBOR::Forward::VolGen::SABR<Process::RandomVec::StdBrownAntithetic>
        lVolGen( inInitVol, inExponent, inVolVol, inCorrSV, inNPath, lTerms,
                 lTenor );
    StepCalculator_ lStep( lTerms, lTenor, inCorr, lVolGen );
    auto lFR =
        LIBOR::Forward::Factory( inNPath, lTerms, lTenor, inInitFR, lStep )
            .template createForwardRates<
                Process::RandomVec::StdBrownAntithetic>();

    Math::Vec lImpVolByCaplet   = Math::makeVec( inIndTenor.size() - 1 );
    Math::Vec lImpVolByFloorlet = Math::makeVec( inIndTenor.size() - 1 );
    std::cout << "Price of Caplet   : " << lFR.calcCaplet( 0.0, 1 )
              << std::endl;
    std::cout << "Price of Floorlet : " << lFR.calcFloorlet( 0.05, 1 )
              << std::endl;
    for ( std::size_t i = 1; i < inIndTenor.size() - 1; ++i )
    {
        lImpVolByCaplet[i]   = lFR.calcBlackImpVolByCaplet( inInitFR[i], i );
        lImpVolByFloorlet[i] = lFR.calcBlackImpVolByFloorlet( inInitFR[i], i );
        for ( double strike = 0.045; strike < 0.065; strike += 0.001 )
        {
            std::cout << lFR.calcBlackImpVolByCaplet( strike, i ) << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "implied volatility by caplet   : ",
        Math::print( lImpVolByCaplet );
    std::cout << "implied volatility by floorlet : ",
        Math::print( lImpVolByFloorlet );
    return 0.0;
}

TEST( ShortRateConstantTest, PriceZCB )
{
    std::vector<std::size_t> lIndTenor{ 0, 10, 20 };
    Math::Vec lInitFR{ 0.05, 0.05 };
    Math::Vec lInitVol{ 0.05, 0.05 };
    Math::Mat lCorr = Math::unitMat( 2, 1.0 );
    Math::Vec lCorrSV{ 0.1, 0.1 };
    double lExponent   = 0.8;
    double lVolVol     = 0.15;
    std::size_t lNPath = 2000000;
    double lMaturity   = 2.0;

    EXPECT_NEAR( 0.0,
                 testImpVolSABR<LIBOR::Forward::StepCalc::NormalTerminalMeas>(
                     lNPath, lMaturity, lIndTenor, lCorr, lInitFR, lInitVol,
                     lExponent, lVolVol, lCorrSV ),
                 0.001 );
}
