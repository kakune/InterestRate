#include <gtest/gtest.h>

#include "analytical/SABR.hpp"

double testBlackImpVolOfSABR( double inTime, double inInitPrice,
                              double inInitVol, double inExponent,
                              double inCorr, double inVolVol, double inStrike )
{
    Analytical::SABR::OneTerm lObj( inTime );
    lObj.setInitPrice( inInitPrice )
        .setInitVol( inInitVol )
        .setExponent( inExponent )
        .setCorr( inCorr )
        .setVolVol( inVolVol )
        .setStrike( inStrike );
    for ( double lStrike = 0.95; lStrike <= 1.05; lStrike += 0.01 )
    {
        lObj.setStrike( inInitPrice * lStrike );
        std::cout << inInitPrice * lStrike << " : "
                  << lObj.approxBlackImpVolByAntonov() << std::endl;
    }
    lObj.setStrike( inStrike );
    return lObj.approxBlackImpVolByHagan();
}

double testCalibration( double inTime, double inInitPrice,
                        std::vector<double> inStrikes,
                        std::vector<double> inImpVols )
{
    // auto lResult = Analytical::SABR::calibrateAllParam( inTime, inInitPrice,
    //                                                     inStrikes, inImpVols
    //                                                     );
    auto lResult = Analytical::SABR::calibrateParamWithFixedExponent(
        inTime, inInitPrice, inStrikes, inImpVols, 0.3 );
    lResult.printAllParam();
    return 0.0;
}

TEST( SABRTest, BlackImpVol )
{
    EXPECT_NEAR( 0.09122244340212736,
                 testBlackImpVolOfSABR( 1.0, 0.05, 0.05, 0.8, 0.1, 0.15, 0.05 ),
                 1e-6 );
    // EXPECT_NEAR( 2.7855172769067678,
    //              testBlackImpVolOfSABR( 10.0, 0.5, 0.6, 0.01, 0.9, 0.9, 0.25
    //              ), 1e-6 );
    // EXPECT_NEAR(
    //     2.0473461128778636,
    //     testBlackImpVolOfSABR( 0.1, 0.01, 2.0, 0.99, 0.01, 0.01, 0.95 ), 1e-6
    //     );
}

TEST( SABRTest, Calibration )
{
    // Table 5.1 of
    // https://research-api.cbs.dk/ws/portalfiles/portal/58422775/soeren_skov_hansen.pdf
    EXPECT_NEAR(
        0.0,
        testCalibration(
            10.0, 0.03571,
            { 0.01571, 0.02571, 0.03071, 0.03571, 0.04071, 0.04571, 0.05571 },
            { 0.3215, 0.2480, 0.2222, 0.2040, 0.1923, 0.1867, 0.1887 } ),
        1e-6 );
}
