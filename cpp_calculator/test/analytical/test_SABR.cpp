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
    for ( double lStrike = 0.75; lStrike <= 1.25; lStrike += 0.02 )
    {
        lObj.setStrike( inInitPrice * lStrike );
        std::cout << lObj.approxBlackImpVolByHagan() << std::endl;
    }
    lObj.setStrike( inStrike );
    return lObj.approxBlackImpVolByHagan();
}

TEST( SABRTest, BlackImpVol )
{
    EXPECT_NEAR( 0.09122244340212736,
                 testBlackImpVolOfSABR( 1.0, 0.05, 0.05, 0.8, 0.1, 0.15, 0.05 ),
                 1e-6 );
    EXPECT_NEAR( 2.7855172769067678,
                 testBlackImpVolOfSABR( 10.0, 0.5, 0.6, 0.01, 0.9, 0.9, 0.25 ),
                 1e-6 );
    EXPECT_NEAR(
        2.0473461128778636,
        testBlackImpVolOfSABR( 0.1, 0.01, 2.0, 0.99, 0.01, 0.01, 0.95 ), 1e-6 );
}
